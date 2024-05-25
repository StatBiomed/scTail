============================================================
scTail for alternative PAS analysis in single cells
============================================================
|pypi| 

.. |pypi| image:: https://badge.fury.io/py/CamoTSS.svg
       :target: https://pypi.org/project/CamoTSS/

.. image:: https://zenodo.org/badge/497821671.svg
      :target: https://zenodo.org/badge/latestdoi/497821671


Note
============
Hi there, my github account did not notify me when there are issue. 
So if you are in a hurry, you can email me. ruiyan@connect.hku.hk.
I check email every day.  



Installation
============

You can install from this GitHub repository for latest (often development) 
version by following command line

.. code-block:: bash

  pip install -U git+https://github.com/StatBiomed/scTail

In either case, add ``--user`` if you don't have the write permission for your 
Python environment.


Quick start
===========

Download test file
===================

You can download test file from figshare_.

.. _figshare: https://figshare.com/articles/dataset/scTail_supplementary_data/25902508

Here, you can download test data and also gene and PAS expression profiles for three dataset: human intestinal, mouse forelimb and ESCC.
  
Run scTail
=============

Here are three steps in scTail : **scTail-callPeak**, **scTail-peakMerge** and **scTail-count**.

We set these three steps to speed up when running some large file (file size > 30G).

Please check your reads1 (the one that contains cellbarcode and UMI) at first before you run scTail to make sure the length of it more than 100bp. In the most situations, it is perfect that length of reads 1 is 150bp or 151bp.

scTail only support two species: mouse and human. Because classifier embedded in it only trains with sequence of mouse and human.

When you get fastq file, you should follow this instruction_ to run scTail step by step. 

.. _instruction: 



Differential APA usage detecting
=================================

To identify differential alternative PAS usage, BRIE2 (Huang & Sanguinetti, 2021) is recommend to be used. 

Here, we provide an example exploiting BRIE2 to detect differential PAS usage. 

You can check it in our manual_.

.. _manual: https://camotss.readthedocs.io/en/latest/runBRIE.html  


Detailed Manual
================

The full manual is here_, including:

`Preprocess`_

`Run CamoTSS`_

`Detect alternative TSS/CTSS`_

.. _here: https://camotss.readthedocs.io/en/latest/index.html

.. _Preprocess: https://camotss.readthedocs.io/en/latest/preprocess.html

.. _Run CamoTSS: https://camotss.readthedocs.io/en/latest/run_CamoTSS.html

.. _Detect alternative TSS/CTSS: https://camotss.readthedocs.io/en/latest/runBRIE.html



Reference
===========














