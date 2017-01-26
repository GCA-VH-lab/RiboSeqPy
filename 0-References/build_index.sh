# build indexes for bowtie2 and hisat2
# Python scripts are part of hisat tutorial vt. Perta M et al, Nature Protocol (2016)
extract_splice_sites.py  Genome.gtf  > Genome.ss
extract_exons.py Genome.gtf  > Genome.exons
bowtie-build ncRNA.fa Indexes/ncRNA
bowtie2-build ncRNA.fa Indexes/ncRNA
hisat2-build ncRNA.fa Indexes/ncRNA
bowtie2-build Genome.fa Indexes/Genome
hisat2-build --ss Genome.ss  --exon Genome.exons Genome.fa  Indexes/Genome

# tabix index
# tabix is distributed with pysam 
grep ^# Genome.gtf; grep -v ^# Genome.gtf | sort -k1,1 -k4,4n | bgzip > genome.gtf.gz
tabix -p gff genome.gtf.gz
