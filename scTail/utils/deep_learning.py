import pickle
import pandas as pd
from scipy import stats
from pyfaidx import Fasta
import pyranges as pr
from tqdm import tqdm
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
from kipoiseq.transforms.functional import one_hot_dna
import numpy as np
import torch
import torch.nn as nn





class scDataset(Dataset):
    def __init__(self,ref_fastq,bed,chromosizedf):
        
        #read bed file
        genes=Fasta(ref_fastq)
        bedfile=pd.read_csv(bed,delimiter='\t')
        #print(len(bedfile))


        bedfile['start']=bedfile['PAS']-100
        bedfile['end']=bedfile['PAS']+100

        #filter to make sure get interger sequence
        bedfile=bedfile[bedfile['start']>1]

        chromosizedf.reset_index(inplace=True)
        chromosizedf.columns=['Chromosome','length']
        bedfile=bedfile.merge(chromosizedf,on='Chromosome')
        bedfile=bedfile[bedfile['PAS']+100 < bedfile['length']]
        bedfile.reset_index(inplace=True,drop=True)
        #print(len(bedfile))



        #extract sequence and do one-hot encoding 
        onehotls=[]
        for chrom, start, end, strand in tqdm(bedfile[['Chromosome','start','end','Strand']].values,desc='Loading Loci'):
            if strand=='+':
                seq_=genes.get_seq(chrom,start+1,end).seq
            else:
                seq_=genes.get_seq(chrom,start+1,end,rc=True).seq

            onehotnp=one_hot_dna(seq_, ['A','C','T','G'])
            onehotls.append(onehotnp)

            
            
        
        # get corresponding X value 
        self.combine_array=np.stack(onehotls)
        self.gene_id=bedfile['cluster_id']

                
    def __len__(self):
        return len(self.combine_array)
    
    
    def __getitem__(self,index):

        seq=np.transpose(self.combine_array[index]).astype(np.float32)
        #print(seq)
        PAS_name=self.gene_id[index]
        #print(PAS_name)
        
        sample = {"seq": seq,'PAS_name':PAS_name}
        return sample



# Model
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv = nn.Sequential(
            nn.Conv1d(4, 128, 8),
             nn.ReLU(),
            nn.Conv1d(128,64,4),
            nn.ReLU(),
            nn.Conv1d(64,32,2),
             nn.BatchNorm1d(32),
              nn.MaxPool2d(2),
            nn.Dropout(0.4),
             nn.Flatten(),
            nn.Linear(1504, 32),
            nn.ReLU(),
        )
        self.fc = nn.Sequential(
            nn.Linear(32, 2),
            nn.Sigmoid()
        )
        
    def forward(self, x):
        x = self.conv(x)
        x = self.fc(x)
        return x





def evaluate(device,net,dataloader):
    net.eval()
    pas_name_ls=[]
    pred_labels=[]
    y_scorels=[]
    with torch.no_grad():
        for data in dataloader:
            seq=data['seq'].to(device)
            output=net(seq)
            #print(output)
            preds = torch.argmax(output, 1)
            pred_labels.append(preds.cpu().data.numpy())

            pas_name=data['PAS_name']
            pas_name_ls.append(pas_name[0])
            y_scorels.append(output[:,1].cpu().data.numpy())
            
    return pas_name_ls,pred_labels,y_scorels



def test(testdataloader,device,net,pretrained_model_path,save_file):
    print(device)

    #os.environ["CUDA_VISIBLE_DEVICES"] = str(device)
    

    check_point=torch.load(pretrained_model_path,map_location=torch.device(device))
    net.load_state_dict(check_point['state_dict'])

    #get predicted label
    pas_name_ls,pred_labels,y_scorels=evaluate(device,net,testdataloader)


    # save prediction
    outputdf=pd.DataFrame({'PAS':pas_name_ls,'y_pred':pred_labels,'y_score':y_scorels})
    outputdf.to_csv(save_file,sep='\t',index=None)
    positivedf=outputdf[outputdf['y_pred'].apply(lambda x: x[0]==1)]
    
    return positivedf



