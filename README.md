# RiboSeqPy

> Transferred from tmargus [github page](https://github.com/tmargus/RiboSeqPy).  
> Changed to organisation.  

Collection of routines for processing **Ribo-Seq** data joined into Python pipeline. Our main aim is to analyze changes of translational dynamic between two different conditions in yeast. We estimate differences in codon occupancy between two conditions and looking for context what is likely related to higher/lower codon occupancy.

The first part of pipeline  starts with _fastq_ preprocessing continues with aligning reads to genome, mapping ribosome positions (uncorrected) and ends with producing metagene plots around start and stop codons. 

Second part of pipeline corrects mapped RPF positions according given offsets in the `readlength_offsets_5-End.txt` to 5' position of P-Site codon: P-Site assignment. It calculates `codon relative rpm` and `codon relative fold difference` (FD) and adds sequence information like codons in E,P,A-Site; nucleic acid sequence and it's translation extending from A-Site to tunnel (by default 11 amino acids). __FD__ calculation assumes you have Ribo-Seq data for two conditions (A condition; B wild type). Important limitation is that __FD__ calculation can't handle multi-exon genes. This is how far it goes in the moment. 


### Prerequisite
1) Python:

  It's ok to use system python but have you own local version gives more flexibility. I used Anaconda Python v.3.5 from [Continuum](https://www.continuum.io/downloads). It comes with a bunch of libraries and have a nice  package manager `conda`. Before `conda` is able to install bioinformatic libraries/programs you have add the _bioconda_ channel.
  
``` 
    conda config --add channels bioconda
```
2) [sra-tools](https://github.com/ncbi/sra-tools/wiki/Downloads) 
          
3) [cutadapt](https://cutadapt.readthedocs.io/en/stable/)

```
    conda install cutadapt
```

4) `pigz` - optional if not installed cutadapt falls back to single core mode

5) [HISAT2](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads)

   _version 2.0.5_ or higher.
   Hisat2 trims ends of reads with bad quality by default. That leads to uncorrect mapping of ribosome location. From the version 2.0.5 there is an option to turn this behavior off.
   
6) [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

7) [samtools](https://github.com/samtools/samtools/) 

```
    conda install samtools
```
  
8) [pysam](https://github.com/pysam-developers/pysam)

```
    conda install pysam
```

### Additional data files

  * Genome.fa  - genome sequence in FastA format
  * ncRNA.fa   - non coding RNA in FastA format
  * Genome.gtf - genome annotation in GTF (gff2) format

Other data files are derived based on those three and commands for that are described in the file  `build_index.sh`.
_Saccharomyces cerevisiae_ genome, annotation, ncRNA and indexes are locating in the folder **0-References/**.
Dummy dataset for testing purpose with __1__ milj. reads locates in the folder **1-Raw/**.


### Usage
Pipeline is split in two parts and controlled by parameters in the file **Param.in**. Part 1 runs up to mapping read ends to genome and producing metagene plots (steps 1 -  8 in Param.in). Part 2 creates corrected P-Site assignment and computes codon relative fold differences (FD) of between 2 conditions. In  **Param.in** you can specify steps you want to run, read length range, mapping (5' or 3') etc.. 

```
    python  Pipeline_part_1.py
```

Edit offset file (`readlength_offsets.txt`) according read length and offsets.

```
    python  Pipeline_part_2.py
```

Logfiles are generated for each step and stored in Reports/ folder.

### Limitations
One important limitation is that calculation of __codon relative fold difference__  can't handle multi-exon genes in the moment. This limitation restricts its usage with bacteria and yeast.

#### Hard coded variables
Some variables are hard coded in to python script. 
(a) Adapter sequence for `cutadapt`  is _CTGTAGGCACCATCAAT_. Locates in `Pipeline_part_1.py` function **cutAdapt()**.
(b) Chromosome names (I, II, III, IV, V, VI, VII, VIII, IX, X, XI, XII, XIII, XIV, XV, XVI, Mito) must be compatible with Genome.fa and Genome.gtf.    
Locates in `Pipeline_part_1.py` and `Pipeline_part_2.py` function **yeastChr()**.
(c) GTF file from ensembl, i. e. must contain features `stop_codon` and `start_codon`  

### References
Firs parts of the code and pipeline backbone is based on a code used in  (1) Radhakrishnan, A., et al. Cell (2016)
https://github.com/greenlabjhmi/2016-Cell-Dhh1. Second part for calculating `codon relative fold difference` is Python 3 adaptation similar to (2) Kannan, K., et al. PNAS (2014)  

