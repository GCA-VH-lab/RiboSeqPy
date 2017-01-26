# RiboSeqPy

Collection of routines for processing **Ribo-Seq** data joined into Python pipeline.
Our main aim is to analyze changes of translational dynamic on different conditions in yeast.
The first part of pipeline starts with _fastq_ preprocessing continues with aligning reads
to genome, mapping ribosome positions (uncorrected) and ends with producing metagenomic
plots around start and stop codons. This is how far it goes in the moment. 
New analyses will be added soon.


### Prerequisite
1) Python:

  Its ok to use system python but have you own local version gives more flexibility.
  I used Anaconda Python v.3.5(& 2.7) from [Continuum](https://www.continuum.io/downloads). It comes with a bunch of libraries and have a nice 
  package manager `conda`. Before conda is able to install bioinformatic libraries/programs you have
  add the _bioconda_ channel.
  
    conda config --add channels bioconda
          
2) [cutadapt](https://cutadapt.readthedocs.io/en/stable/)   

    conda install cutadapt
          
3) HISAT2    ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads

   _version 2.0.5_ or higher.
   Hisat2 trims ends of reads with bad quality by default. That leads to uncorrect mapping of ribosome location.
   From the version 2.0.5 there is an option to turn this behavior off.
   
4) [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

5) [samtools](https://github.com/samtools/samtools/) 

    conda install samtools
        
6) [pysam](https://github.com/pysam-developers/pysam)

    conda install pysam

### Additional data files

  * Genome.fa  - genome sequence in FastA format
  * ncRNA.fa   - non coding RNA in FastA format
  * Genome.gtf - genome annotation in GTF (gff2) format

Other data files are derived based on those three and commands for that are described in the file  `build_index.sh`.
_Saccharomyces cerevisiae_ genome, annotation, ncRNA and indexes are locating in the folder **0-References/**.
Dummy dataset for testing purpose with 1 milj reads locates in the folder **1-Raw/**.


### Usage
Pipeline is controlled by parameters in the file **Param.in**. You can specify steps you want to run,
readlength range, mapping (5' or 3'). 

    python  Pipeline_iv.py

Logfiles are generated for each step and stored in Reports/ folder.


### References
Firs parts of the code and pipeline backbone is based on a code used in Radhakrishnan, A., et al. Cell (2016)
https://github.com/greenlabjhmi/2016-Cell-Dhh1
