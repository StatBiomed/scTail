from ..version import __version__

def main():
    print("Welcome to scTail v%s! Command lines available:\n " %(__version__))
    print("scTail-callPeak\n Get PAS cluster from bam file of each samples")
    print("scTail-peakMerge\n Merge PAS clusters from different samples")
    print("scTail-count \n Assign reads2 to PAS clusters")

if __name__ == "__main__":
    main()