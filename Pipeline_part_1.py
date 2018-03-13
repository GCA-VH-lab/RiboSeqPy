# For Python 2 comment these 2 lines out
#from __future__ import division
#from __future__ import print_function
#
# run code:
#       python Pipeline_part_1.py
#
import os
import time
import gzip
import shutil
import pysam
import subprocess as sp
import numpy as np # optional,  used once
import numpy.random as rn
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')      # load backend - server safe
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from collections import defaultdict

import warnings
import tables

# to supress NaturalNameWarning rawAssignment()
# pd.DataFrame column names are '25', '26' ... etc and saving to h5 will cause warning
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

#  Codes backbone from Radhakrishnan, A., et al. Cell (2016)
#!  Changes
#! adapted for Python 3
#! Mapping part is rewritten -> using pysam
#! Annotation (GTF format) is managed by tabix (comes with pysam)
#! Metagenomic tables and plots were added
#! SamtoBam(iX) *.sam -> sorted & indexed *.bam
#
# This code requires the following programs to be installed on your computer
#
# 1) sra-tools (https://github.com/ncbi/sra-tools/wiki/Downloads)
# 2) cutadapt  (https://cutadapt.readthedocs.io/en/stable/)   #	v 1.15 parallel version    % conda install cutadapt
# 3) hisat2    (ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads)   version 2.0.5 or higher
# 4) bowtie2   (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
# 5) samtools  (https://github.com/samtools/samtools/)        #	via anaconda   % conda install samtools
# 6) pysam     (https://github.com/pysam-developers/pysam)    # via anaconda   % conda install pysam
# 7) pigz      (https://github.com/madler/pigz)               # to work in pair with parallel cutadapt
# 8)
# qualityFilter applies filtering for readlength reading rage from Param.in
#

def cleanFile(File, Condition):

    # gZip a file and delete the un-gZipped version!
    if Condition == "gzip":
        with open(File, 'rb') as FIn, gzip.open(File + ".gz", 'wb') as FOut:
            shutil.copyfileobj(FIn, FOut)

        sp.Popen(["rm", File])

    if Condition == "bgzip":
        BgZip = ["bgzip", "-@", Params['cpu'], File]
        BgZipIt = sp.Popen(BgZip, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
        BgZipIt.wait()

    if Condition == "rm":
        sp.Popen(["rm", File])


def makeDirectory(Path):

# Check if a folder named exists at Path. If not, create it!
    
    Split = Path.split("/")
    
    if len(Split) == 1:
        List = sp.Popen(["ls"], stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True) 
    else:
        List = sp.Popen(["ls", "/".join(Split[:-1])], stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)

    if Path in List.communicate()[0].split("\n"):
        pass
    else:
        Make = sp.Popen(["mkdir", Path])
        Make.wait()


def SamToBam(iX):
    
# sort and convert SAM -> BAM + index

    SamFile     = "5-Aligned/" + iX + ".sam"
    BamFile     = "5-Aligned/" + iX + ".bam"
    
    SamTools    = ["samtools", "sort", "-m 4G", "-@ 4", "-o", BamFile, SamFile]
    SortConvert = sp.Popen(SamTools, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True) 
    SortConvert.wait()
    
    SamTools    = ["samtools", "index",  BamFile]
    IndexBam    = sp.Popen(SamTools, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
    IndexBam.wait()
    
    cleanFile(SamFile, "rm")


def parseParams(Path):

# Open the parameter file and read in parameters and files to operate on

    File      = open(Path)
    SRAList   = []
    Names  = []
    ParamDict = {}

    for iL in File:
        if iL[0] == "#":
            Split = iL[2:-1].split(":")
            if len(Split) > 1:
                ParamDict[Split[0].strip()] = Split[1].strip()
        else:
            if len(iL) != 1:
                Split = iL[:-1].split("\t")
                SRAList.append(Split[0])
                Names.append(Split[1])

    return ParamDict, SRAList, Names


def downloadData(SRAList, NameList, Params):
    # Check to see if all the files
    makeDirectory("1-Raw")

    for iX in range(len(SRAList)):
        Dump  = sp.Popen(["fastq-dump", SRAList[iX]], stdout=sp.PIPE, stderr=sp.PIPE)
        Dump.wait()

        Move  = sp.Popen(["mv", SRAList[iX] + ".fastq", "1-Raw/" + NameList[iX] + ".fastq"])
        Move.wait()
        cleanFile("1-Raw/" + NameList[iX] + ".fastq", Params["Clean"]) # original
        #cleanFile("1-Raw/" + NameList[iX] + ".fastq", "bgzip")


def update_df(df, Chr, strand):
    df.fillna(0, inplace=True)
    df["sum"] = df.sum(axis=1)
    columns = list(df.columns)
    columns = ["Chr", "Position", "Strand"] + columns
    df["Chr"] = Chr
    df["Strand"] = strand
    df["Position"] = df.index

    return df[columns]


def not_enough_data(df, Threshold=12):
    '''True when Not enough data
    Assumes numeric columns only
    Project: Yeast-eEF3
    '''
    return True if df.sum().sum() < Threshold else False


def raw_metag_threshold_to_rpm(BamName, Threshold):
    """
    Converts raw MetagThresold to rpm MetagThreshold

    :param BamName:    BAM file used for caclulating normalization factor
    :param Threshold:  raw threshold
    :return: normalized threshold to fit with rpm normlized data
    """
    bamfile = pysam.AlignmentFile(BamName, "rb")  # open BAM file
    c = 0
    for read in bamfile.fetch():
        if read.get_tag("NH") == 1: #mapped once
            c+=1

    return Threshold/(c/10**6)


def df_framing(df1, index, columns, strand="+"):
    """ returns df what contains values for all positions in the given range.
    Original df is condensed, i. e. positions with values < 0 not included

    :param df1:     condensed df, i. e.  don't contain rows with 0 in index
    :param index:   range of genome positions
    :param columns: list of read length + 'sum'
    :param strand:  default "+"  alt. "-"
    :return: merged df
    """
    # create df2
    df2 = pd.DataFrame(0, index=index, columns=columns)
    df1 = df1.add(df2, fill_value=0, axis=1)
    if strand == "+":
        df1.reset_index(inplace=True)  # reset index
        return df1[columns]
    elif strand == "-":
        df1 = df1[::-1]  # reverts table
        df1.reset_index(inplace=True)  # reset index
        return df1[columns]
    else:
        # error
        print("ERROR! Expext '+'/'-' but found {} for strand".format(strand))


def dfTrimmiX3(df, Span, iX):
    """Truncates Data Frame to fit in figure 3pr"""
    #todo: rewrite it similarly to dfTrimmiX5()
    x = 3 * 6
    if (Span == 60) & (iX == "Start"):
        df = df[30:90 + 2 * x]
        df.reset_index(inplace=True)
        del df['index']
    elif (Span == 60) & (iX == "Stop"):
        df = df[30:90 + 2 * x]
        df.reset_index(inplace=True)
    else:
        pass

    return df


def dfTrimmiX5(df, Span, iX, inside_gene=33, outside_gene=18 ):
    """Truncates Data Frame to fit in figure 5pr """
    if (inside_gene > Span) | (outside_gene > Span):
        print("Given parameters inside- or outside gene are bigger than Span!\nQuering out of range data!")
        return df

    if iX == "Start":
        return df.loc[-outside_gene:inside_gene, ]
    elif iX == "Stop":
        return df.loc[-inside_gene:outside_gene, ]
    else:
        print("Table is not modified. Mapping is unkown!")
        return df


def colorsCheck(dic, key):
    # import numpy.random as rn
    '''Generates random color for missing key'''
    if key not in dic.keys():
        a = rn.rand(3, 1)
        l = [a[i][0] for i in list(range(3))]  # decomposing array
        dic[key] = tuple(l)  # list 2 tuple

    return dic


def restructurate_hd5(infile, outfile, close_outfile=True):
    """ infile.h5 keys - "/For_raw", "/Rev_raw", ...
    outfile.h2  keys - "/For_raw/I", "/For_raw/II", ... etc
    "Position" is set to index

    :param infile:
    :param outfile:
    :return: reindexed 2 level hdf
    """
    # open inp_HDF
    inp__h5 = pd.HDFStore(infile, "r")
    outp_h5 = pd.HDFStore(outfile, complevel=5, complib="zlib", mode="w")
    # open out_HDF

    # for each I level table For_raw, Rev_raw, ...
    for key in inp__h5.keys():
        # for each chromosome
        df = inp__h5[key]
        for Ch in df['Chr'].unique():
            df_ch = df[df.Chr == Ch].copy()

            # set Position to index
            df_ch.set_index('Position', inplace=True)

            # save df under II level key what is now chromosome
            h5_key = key + "/" + Ch
            outp_h5.put(h5_key, df_ch)

    inp__h5.close()

    # todo: close file in function and open it later
    # leaving it open here will couse problems later, because it is not

    if close_outfile == True:
        outp_h5.close()
    else:
        return outp_h5


def yeastChr():
    # Ordered yeast Chr list short names from ensembl
    return ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','Mito']

# generalized do   "Bioinformatics Programming Using Python by Mitchell L Model"
def do(collection, fn):
    '''  Generalized do function
    '''
    for item in collection:
        fn(item)

# redefine print   "Bioinformatics Programming Using Python by Mitchell L Model"
def print_collection(collection):
    ''' Generalized print_collection function
    '''
    do(collection, print)

def print_params(Params):
    print_collection(("{:15s}- {}".format(k, Params[k]) for k in sorted(Params)))


def cutAdapt(SRAList, Names, Params):

    makeDirectory("2-Trimmed")
    makeDirectory("2-Trimmed/Reports")

    for iX in Names:
        #eEF3
        CutAdapt     = ["cutadapt","--discard-untrimmed", "-j", Params['cpu'], "-a","CTGTAGGCACCATCAAT","-o","2-Trimmed/","1-Raw/"]
        CutAdapt[7] += iX + ".fastq.gz"
        CutAdapt[8] += iX + ".fastq.gz"
        #New1
        #CutAdapt = ["cutadapt", "--discard-untrimmed", "-j", Params['cpu'], "-u", "3","-a", "NNNNCTGTAGGCACCATCAAT","-o", "2-Trimmed/", "1-Raw/"]
        #CutAdapt[9] += iX + ".fastq.gz"
        #CutAdapt[10] += iX + ".fastq.gz"
        report = "cutadapt command line:  {}".format(CutAdapt); print(report)
        Trim         = sp.Popen(CutAdapt, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
        Trim.wait()

        FileOut = open("2-Trimmed/Reports/" + iX + ".txt","w")
        FileOut.write(Trim.communicate()[0])
        FileOut.close()
        

def qualityFilter(SRAList, Names, Params):
    
    PHREDDict = {
        "!": 9.999999e-01, "\"": 7.943282e-01, "#": 6.309573e-01, "$": 5.011872e-01, "%": 3.981072e-01,
        "&": 3.162278e-01, "\'": 2.511886e-01, "(": 1.995262e-01, ")": 1.584893e-01, "*": 1.258925e-01,
        "+": 1.000000e-01, ",": 7.943282e-02, "-": 6.309573e-02, ".": 5.011872e-02, "/": 3.981072e-02,
        "0": 3.162278e-02, "1": 2.511886e-02, "2": 1.995262e-02, "3": 1.584893e-02, "4": 1.258925e-02,
        "5": 1.000000e-02, "6": 7.943282e-03, "7": 6.309573e-03, "8": 5.011872e-03, "9": 3.981072e-03,
        ":": 3.162278e-03, ";": 2.511886e-03, "<": 1.995262e-03, "=": 1.584893e-03, ">": 1.258925e-03,
        "?": 1.000000e-03, "@": 7.943282e-04, "A": 6.309573e-04, "B": 5.011872e-04, "C": 3.981072e-04,
        "D": 3.162278e-04, "E": 2.511886e-04, "F": 1.995262e-04, "G": 1.584893e-04, "H": 1.258925e-04,
        "I": 1.000000e-04, "J": 7.943282e-05
    }

    makeDirectory("3-Filtered")
    makeDirectory("3-Filtered/Reports")
    LogFileName = "3-Filtered/Reports/" + "Quality_filtering" + "_iv_log.txt"
    LOG_FILE = open(LogFileName, "wt")

    for iX in Names:
        low_qual     = 0 
        short_reads  = 0
        long_reads   = 0
        Quality      = float(Params["Quality"])
        read_min_len = int(Params["ReadLenMiN"]) #- fshift
        read_max_len = int(Params["ReadLenMaX"]) #+ fshift
        # done: crashes when 2-Trimmed/Reports/*.txt does not exists - happens when you start already trimmed input
        Length = ""
        cutadapt_report = "2-Trimmed/Reports/" + iX + ".txt"

        if not os.path.exists(cutadapt_report):
            os.system("gzcat  2-Trimmed/" + iX + ".fastq.gz |wc -l > tmp.t")
            lines_in_fastq = int(open('tmp.t', 'r').read()[:-1])
            Length = int(lines_in_fastq / 4)
            print("{} adapters removed elsewhere! \n".format(iX))
        else:
            File = open(cutadapt_report)
            # todo: cutadapt version inconsistency - new version have 1 more line _- range(0,9) instead (0,8)
            Burn = [File.readline() for Idx in range(0, 9)]
            Length = int(Burn[-1][:-1].split(" ")[-1].replace(",", ""))
            File.close()

        report = " {:18} : {:>12,}".format(iX, Length); print(report)
        LOG_FILE.write(report + "\n")

        File    = gzip.open("2-Trimmed/" + iX + ".fastq.gz", "rt")
        FileOut = open("3-Filtered/" + iX + ".fastq", "w")

        for iN in range(0,Length):
            Identifier  = File.readline().rstrip("\n")
            Sequence    = File.readline().rstrip("\n")
            QIdentifier = File.readline().rstrip("\n")
            PHRED       = File.readline().rstrip("\n")
            Score       = 1.0
            Len         = len(PHRED)
            
            if  Len < read_min_len:
                short_reads += 1

            elif Len > read_max_len:
                long_reads  += 1

            else:
                for IdxL in range(0,Len):
                    Score   = Score*(1 - PHREDDict[PHRED[IdxL]])

                if (Score > Quality):
                    FileOut.write(Identifier + "\n" + Sequence + "\n" + QIdentifier + "\n" + PHRED + "\n")
                else:
                    low_qual +=1

        #SumPhred = (-10 * np.log10(1 - Quality), low_qual) # Convert Quality to Phred sum

        report  = " Reads len < {:>3} : {:>12,}\n".format(read_min_len, short_reads)
        report += " Reads len > {:>3} : {:>12,}\n".format(read_max_len, long_reads)
        report += " Quality   < {:>3} : {:>12,}\n".format(Quality, low_qual)
        report += " Reads left      : {:>12,}".format(Length - (short_reads + low_qual))

        print(report, "\n"); LOG_FILE.write(report + "\n\n")

        File.close()
        FileOut.close()

    LOG_FILE.close()


def ncRNASubtract(SRAList, Names, Params):
    
    makeDirectory("4-Subtracted")
    makeDirectory("4-Subtracted/SAM")
    makeDirectory("4-Subtracted/Reports")

    for iX in Names:
        ncRNA    = "0-References/Indexes/ncRNA"
        Input    = "3-Filtered/" + iX + ".fastq"
        Output   = "4-Subtracted/SAM/" + iX + ".sam"
        Unmapped = "4-Subtracted/" + iX + ".fastq"
        
        Bowtie2  = ["bowtie2", "--no-unal", "-p", Params['cpu'], "--un", Unmapped, "-x", ncRNA, "-U", Input,"-S", Output]
        print(Bowtie2)
        Subtract = sp.Popen(Bowtie2, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
        Subtract.wait()
        
        FileOut  = open("4-Subtracted/Reports/" + iX + ".txt","w")
        FileOut.write(Subtract.communicate()[1])
        FileOut.close()
        print("Cleaning: {} \t{}".format("3-Filtered/" + iX + ".fastq",Params["Clean"]))
        cleanFile("3-Filtered/" + iX + ".fastq", Params["Clean"])
        cleanFile("4-Subtracted/SAM/" + iX + ".sam", "rm")


def genomeAlign(SRAList, Names, Params):

    makeDirectory("5-Aligned")
    makeDirectory("5-Aligned/Reports")

    for iX in Names:
        Genome  = "0-References/Indexes/Genome"
        Input   = "4-Subtracted/" + iX + ".fastq"
        Output  = "5-Aligned/" + iX + ".sam"
        NotAli  = "5-Aligned/" + iX + "_unaligned.fastq"  # added! unalinged reads output
        
        Hisat2  = ["hisat2", "--no-unal", "-p", Params['cpu'], "--no-softclip", "--dta", "-x", Genome, "-U", Input, "-S", Output]
        # IF reads not aligned to genome also needed ncomment line bleow
        #Hisat2 = ["hisat2", "-p 6", "--no-softclip", "--dta", "--un", NotAli, "-x", Genome, "-U", Input, "-S", Output]
        print(Hisat2)

        Align   = sp.Popen(Hisat2, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)

        FileOut = open("5-Aligned/Reports/" + iX + ".txt", "w")
        FileOut.write(Align.communicate()[1])
        FileOut.close()

        cleanFile("4-Subtracted/" + iX + ".fastq", Params["Clean"])
        print("Sam2Bam {}".format(iX))
        SamToBam(iX)


def rawAssignment(SRAList, Names, Params):
    # changes: 04.Jul :: restucturate h5 file - creates h5 files with keys like "For_rpm/I" and sets "Position" to index
    makeDirectory("6-AssignRaw")
    makeDirectory("6-AssignRaw/Reports")
    # include_mapped_twice influence ho normalisation factor is computed. vt blow
    include_mapped_twice = Params['MappedTwice']  # includes reads mapped twice NH:i:2
    save_csv             = False  # save output to tab delim csv. in addition to hdf5

    rlmin   = int(Params["ReadLenMiN"])
    rlmax   = int(Params["ReadLenMaX"])
    Mapping = Params["Mapping"]  # Mapping 5 or 3 prime end
    rlrange = str(rlmin) + "-" + str(rlmax) # read length range 4 filename

    for iN in Names:

        BamName     = "5-Aligned/" + iN + ".bam"          # sorted and indexed BAM
        bamfile     = pysam.AlignmentFile(BamName, "rb")  # open BAM filee
        outfile_for = "6-AssignRaw/" + iN + "_" + Mapping + "-End_" + rlrange + "_raw_iv" + "_For.txt"
        outfile_rev = "6-AssignRaw/" + iN + "_" + Mapping + "-End_" + rlrange + "_raw_iv" + "_Rev.txt"
        outfile_hdf = "6-AssignRaw/" + iN + "_" + Mapping + "-End_" + rlrange + "_iv" + ".h5"
        outf_idx_hdf= "6-AssignRaw/" + iN + "_" + Mapping + "-End_" + rlrange + "_idx_iv" + ".h5"
        LogFileName = "6-AssignRaw/Reports/" + iN + "_" + Mapping + "-End_" + rlrange + "_iv_log.txt"
        LOG_FILE    = open(LogFileName, "wt")
        # counters for log
        c2_twice = c_once = total_no = 0
        # empty dataframe for collecting data
        df_for_sum = pd.DataFrame()
        df_rev_sum = pd.DataFrame()
        # Process Log
        report = "\nBamFile: {}\nrlmin:   {}\nrlmax:   {}\nName:    {}\nMapping: {}".format(BamName, rlmin, rlmax, iN,
                                                                                            Mapping)
        LOG_FILE.write(report + "\n"); print(report, "\n")

        # humanChr() gives an ordered list
        for ref in yeastChr():
            c1  = 0
            c2  = 0
            reads_mapped_ref   = 0
            ref_total_read_count = 0

            defF = defaultdict(list)  # DefaultDict  For
            defR = defaultdict(list)  # DefaultDict  Rev
            ForDict = {}  # Collecting  data For
            RevDict = {}  # Collecting  data Rev

            for read in bamfile.fetch(ref):
                ref_total_read_count += 1
                # collect no of readsa
                if (read.get_tag("NH") == 1):
                    c1 += 1
                elif (read.get_tag("NH") == 2):
                    c2 +=1
                else:
                    pass

                if (read.get_tag("NH") == 1):  # NH tag (NH:i:1) tells how many times read are mapped to genome
                    reads_mapped_ref +=1
                    readl = read.query_length  # get read length
                    # Redefining leftmost & rightmost
                    if not read.is_reverse:  # read is Forward
                        beg = read.reference_start    # 5'
                        end = read.reference_end - 1  # 3' correct by -1
                    else:  # read is Reverse
                        beg = read.reference_end - 1  # 5' correct by -1
                        end = read.reference_start    # 3'

                    if Mapping == "5":
                        defR[readl].append(beg) if read.is_reverse else defF[readl].append(beg)
                    if Mapping == "3":
                        defR[readl].append(end) if read.is_reverse else defF[readl].append(end)
                # to include those mapped twice
                if (read.get_tag("NH") == 2) & (include_mapped_twice == "Yes"):
                    reads_mapped_ref += 1
                    readl = read.query_length  # get read length
                    # Redefining leftmost & rightmost
                    if not read.is_reverse:  # read is Forward
                        beg = read.reference_start  # 5'
                        end = read.reference_end - 1  # 3' correct by -1
                    else:  # read is Reverse
                        beg = read.reference_end - 1  # 5' correct by -1
                        end = read.reference_start  # 3'

                    if Mapping == "5":
                        defR[readl].append(beg) if read.is_reverse else defF[readl].append(beg)
                    if Mapping == "3":
                        defR[readl].append(end) if read.is_reverse else defF[readl].append(end)

            dummy = [0]
            for rlen in range(rlmin, rlmax + 1):
                ForDict[rlen] = Counter(defF.get(rlen, dummy))  # .get() method if rlen
                RevDict[rlen] = Counter(defR.get(rlen, dummy))  # if don't exist use dummy

            df_for = update_df(pd.DataFrame(ForDict), Chr=ref, strand="+")
            df_rev = update_df(pd.DataFrame(RevDict), Chr=ref, strand="-")

            df_for_sum = pd.concat([df_for_sum, df_for], ignore_index=True) # collect summary table
            df_rev_sum = pd.concat([df_rev_sum, df_rev], ignore_index=True) # collect summary table


            # Log_File pr Chr
            report = "{:<5s}\t{:>10,d} reads".format(ref, reads_mapped_ref)
            LOG_FILE.write(report + "\n"); print(report)
            # Reset/collect counter data
            c_once   += c1
            c2_twice += c2 # mapped twice
            total_no += ref_total_read_count
            reads_mapped_ref = ref_total_read_count = 0

        # Per Name !!!
        # convert int column names to str
        df_for_sum.rename(columns={i: str(i) for i in range(rlmin, rlmax + 1)}, inplace=True)  # num col_names to str
        df_rev_sum.rename(columns={i: str(i) for i in range(rlmin, rlmax + 1)}, inplace=True)  # num col_names to str

        ## Log Report summary
        report = "\nTotal No of reads {:>11,} mapped to genome\n".format(total_no)
        report += "Number   of reads {:>11,d} mapped once to genome\n".format(c_once)
        report += "Number   of reads {:>11,d} mapped twice to genome reports # will be added if MappedTwice == True \n".format(c2_twice)
        report += "Number   of reads {:>11,d} mapped more than counted already\n".format(total_no - (c_once +c2_twice))
        LOG_FILE.write(report); print(report)
        ##>>
        report = "\nOutput tables are stored:"; LOG_FILE.write(report+ "\n");print(report)

        if save_csv == True:
            df_for_sum.to_csv(outfile_for, sep='\t', header=True, index=True)  # csv table output
            df_rev_sum.to_csv(outfile_rev, sep='\t', header=True, index=True)  # csv table output
            report = "{}\n{}\n".format(outfile_for, outfile_rev)
            LOG_FILE.write(report + "\n")
            print(report)

        report = "{}\tkeys: 'For_raw', 'Rev_raw'".format(outfile_hdf)

        store = pd.HDFStore(outfile_hdf, complevel=5, complib="zlib", mode="w")
        store.put("For_raw", df_for_sum, format="table", data_columns=True)
        store.put("Rev_raw", df_rev_sum, format="table", data_columns=True)
        LOG_FILE.write("\n" + report + "\n"); print(report)
        # Convert to RPM s
        report = "\n Converting raw -> rpm \n"
        LOG_FILE.write(report + "\n"); print(report); report = ""

        ## Convert RAW -> RPM
        # NB! include_mapped_twice = True  influence how RPM is calculated
        normFactor = 0
        if include_mapped_twice == "Yes":
            l = [0 for read in bamfile.fetch() if read.get_tag("NH") <= 2]  # reads mapped once & twice
            normFactor = len(l) / 10 ** 6  # normalisation factor
            report = "Normalization factor {} is computed based on all mapped reads {:,}".format(normFactor, len(l))
        else:
            l = [0 for read in bamfile.fetch() if read.get_tag("NH") == 1]  # reads mapped once
            normFactor = len(l) / 10 ** 6  # normalisation factor
            report = "Normalization factor {} is computed based on reads mapped once {:,}".format(normFactor, len(l))

        LOG_FILE.write(report + "\n"); print(report); report = ""
        col2norm = [str(i) for i in (range(rlmin, rlmax + 1))] + ["sum"]

        for iX in col2norm:  # normalization

            df_for_sum[iX] = df_for_sum[iX] / normFactor
            df_rev_sum[iX] = df_rev_sum[iX] / normFactor
            line = "Normal factor for {} - {:7.4f}".format(iX, normFactor)
            report += line + "\n"
            print(line)

        LOG_FILE.write(report + "\n"); print("")

        if save_csv == True:
            outfile_for, outfile_rev = outfile_for.replace("_raw", "_rpm"), outfile_rev.replace("_raw", "_rpm")
            df_for_sum.to_csv(outfile_for, sep='\t', header=True, index=True)  # csv table output
            df_rev_sum.to_csv(outfile_rev, sep='\t', header=True, index=True)  # csv table output
            report = "{}\n{}\n".format(outfile_for, outfile_rev)
            LOG_FILE.write(report + "\n"); print(report)

        store.put("For_rpm", df_for_sum, format="table", data_columns=True)
        store.put("Rev_rpm", df_rev_sum, format="table", data_columns=True)
        store.close()
        report = "\n{}\tkeys: 'For_rpm', 'Rev_rpm'\n".format(outfile_hdf)
        report += "\n{}\tTime taken thus far: {}".format(iN, time.time() - Start)
        LOG_FILE.write(report + "\n"); print(report)

        # restructurate hdf
        infile  = outfile_hdf
        outfile = outf_idx_hdf
        restructurate_hd5(infile, outfile, close_outfile=True)
        report = "Restructurate hdf\nInfile:   {}\nOutfile    {}".format(infile, outfile)
        LOG_FILE.write(report + "\n"); print(report, "\n")

    LOG_FILE.close()
    bamfile.close()


def metagTables(SRAList, Names, Params):
    makeDirectory("7-MetagTbl")
    makeDirectory("7-MetagTbl/Reports")
    #time.sleep(0.1)
    rlmin = int(Params["ReadLenMiN"])
    rlmax = int(Params["ReadLenMaX"])
    Span = int(Params["MetagSpan"])  # nr of nt before and after 5' position of start/stop codons
    Mapping = Params["Mapping"]  # 5 & 3
    dataNorm = Params["Normalised"]  # "raw" or "rpm"


    columns = [str(i) for i in range(rlmin, rlmax + 1)] + ['sum']
    rlrange = str(rlmin) + "-" + str(rlmax)  # readlength range -> filename
    LogFileName = "7-MetagTbl/Reports/" + "MetagTabl_" + Mapping + "-End_" + rlrange + "_iv_log.txt"
    LOG_FILE = open(LogFileName, "wt")
    #time.sleep(0.1)

    for iN in Names:

        cf1 = cr1 = cf2 = cr2 = 0  # counters
        report  = "\nName: {}".format(iN)
        LOG_FILE.write(report + "\n"); print(report)

        # file names
        fn_body    = iN + "_" + Mapping + "-End_" + rlrange  # filename body
        outf_start = "7-MetagTbl/" + fn_body + "_" + dataNorm + "_Start" + "_iv_Meta_Sum.txt"
        outf_stop  = "7-MetagTbl/" + fn_body + "_" + dataNorm + "_Stop" + "_iv_Meta_Sum.txt"
        infile_h5  = "6-AssignRaw/" +fn_body + "_iv" + ".h5"

        # Empty DataFrames
        meta_start_dff = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)
        meta_start_dfr = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)
        meta_stop_dff  = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)
        meta_stop_dfr  = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)

        report = "Collecting data around genes Start & Stop - Span: {}".format(Span);
        print(report)
        # Annotation
        tabixfile = pysam.TabixFile("0-References/genome.gtf.gz", parser=pysam.asGTF())
        # Adjust Threshold if for RPM if Normalization is "rpm"
        Threshold = int(Params["MetagThreshold"])  # has to be here
        if dataNorm == "rpm":
            BamName = "5-Aligned/" + iN + ".bam"  # sorted and indexed BAM
            Threshold = raw_metag_threshold_to_rpm(BamName, Threshold)
        else:
            pass

        Fkey, Rkey = ("For_rpm", "Rev_rpm") if dataNorm == "rpm" else ("For_raw", "Rev_raw")  # keys for h5
        df_f = pd.read_hdf(infile_h5, Fkey)  # read in Forward str. df from hdf
        df_r = pd.read_hdf(infile_h5, Rkey)  # read in Reverse str. df from hdf

        for Chr in yeastChr():

            dfc_f = df_f[df_f.Chr == Chr].set_index("Position") # get chr subset and index by Position
            dfc_f = dfc_f[columns]                              # selected columns only

            dfc_r = df_r[df_r.Chr == Chr].set_index("Position") # get chr subset and index by Position
            dfc_r = dfc_r[columns]

            for gtf in tabixfile.fetch(reference=Chr):
                # Stop codon in Forw strand
                if (gtf.feature == 'stop_codon') & (gtf.strand == '+'):

                    df = dfc_f.loc[gtf.start - Span:gtf.start + Span][
                        columns].copy()  # subDataFrame  +/-Span around the feature
                    if not_enough_data(df, Threshold=Threshold):  # Check for data /&
                        continue

                    index = range(gtf.start - Span, gtf.start + Span + 1)
                    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
                    meta_stop_dff = meta_stop_dff + df  # sum dataframes
                    cf2 += 1

                # Stop codon in Rev strand
                elif (gtf.feature == 'stop_codon') & (gtf.strand == '-'):
                    df = dfc_r.loc[gtf.end - Span -1:gtf.end + Span -1][columns].copy()  # -1 correction for rev strand
                    if not_enough_data(df, Threshold=Threshold):  # Check for data
                        continue

                    index = range(gtf.end - Span -1, gtf.end + Span)  # -1 correction
                    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
                    meta_stop_dfr = meta_stop_dfr + df  # sum dataframes
                    cr2 += 1

                # Sart codon in Forw strand
                elif (gtf.feature == 'start_codon') & (gtf.strand == '+'):

                    df = dfc_f.loc[gtf.start - Span:gtf.start + Span][
                        columns].copy()  # subDataFrame  +/-Span around the feature
                    if not_enough_data(df, Threshold=Threshold):  # Check for data
                        continue

                    index = range(gtf.start - Span, gtf.start + Span + 1)
                    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
                    meta_start_dff = meta_start_dff + df  # sum dataframes
                    cf1 += 1

                # Start codon in Rev strand
                elif (gtf.feature == 'start_codon') & (gtf.strand == '-'):
                    df = dfc_r.loc[gtf.end - Span -1:gtf.end + Span -1][columns].copy() # -1 correction for rev strand
                    if not_enough_data(df, Threshold=Threshold):  # Check for data
                        continue

                    index = range(gtf.end - Span -1, gtf.end + Span) # -1 correction for rev strand
                    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
                    meta_start_dfr = meta_start_dfr + df  # sum dataframes
                    cr1 += 1

        print("Summing up")
        LOG_FILE.write("Summing up ...\n")
        # summing up
        meta_start_sum = meta_start_dff + meta_start_dfr
        meta_stop_sum = meta_stop_dff + meta_stop_dfr
        # saving to file
        report = "Around START: {}".format(outf_start);
        LOG_FILE.write(report + "\n"); print(report)
        # print("Sum of saved table: {}".format(int(meta_start_sum["sum"].sum())))
        meta_start_sum.to_csv(outf_start, sep='\t', header=True, index=True)
        meta_start_sum['rel_Pos'] = list(range(-Span, Span + 1))
        meta_start_sum.to_csv(outf_start, sep='\t', header=True, index=True)

        report = "Around STOP:  {}".format(outf_stop)
        LOG_FILE.write(report + "\n"); print(report)
        meta_stop_sum['rel_Pos'] = list(range(-Span, Span + 1))
        meta_stop_sum.to_csv(outf_stop, sep='\t', header=True, index=True)

    LOG_FILE.close()


def metagPlots(SRAList, Names, Params):
    # using pandas, matplotlib, seaborn, numpy
    makeDirectory("8-MetagPlot")

    sns.set_style("white")    # seaborn_aesthetic
    sns.set_context("paper")  # seaborn_aesthetic


    Span     = int(Params["MetagSpan"])
    Mapping  = Params["Mapping"]
    dataNorm = Params["Normalised"]  # Mapping 5 or 3 prime end
    rlrange  = Params["ReadLenMiN"] + "-" + Params["ReadLenMaX"]  # readlength range -> filename
    readLen_l= [str(i) for i in range(int(Params["ReadLenMiN"]), int(Params["ReadLenMaX"]) + 1)] + ["sum"]

    # colors for plot
    colors = {'25': 'fuchsia', '26': 'blueviolet', '27': 'darkblue', '28': 'b', '29': 'r',
              '30': 'salmon', '31': 'orange', '32': 'olive', '33': 'g', '34': 'tan',
              '35': 'y', 'sum': 'brown', 'ttl': 'brown'
              }
    a = 0.5  # alphat
    x = 10  # figure width
    y = 1.4 * len(readLen_l)  # figure height

    for iN in Names:
        for iX in ["Start", "Stop"]:
            infile = "7-MetagTbl/" + iN + "_" + Mapping + "-End" + "_" + rlrange + \
                     "_" + dataNorm + "_" + iX + "_iv_Meta_Sum.txt"
            outfig = "8-MetagPlot/" + iN + "_" + Mapping + "-End" + "_" + rlrange + \
                     "_" + dataNorm + "_" + iX + "_iv.png"

            if os.path.isfile(infile):  # infile exits
                print("\n{}".format(outfig))

                fig, axes = plt.subplots(nrows=len(readLen_l), figsize=(x, y))
                fig.suptitle(outfig, y=1.02, fontsize=12)
                df = pd.read_csv(infile, index_col="rel_Pos", sep='\t')
                df.drop("Unnamed: 0", axis=1, inplace=True)

                # Adjust plot for mapping and Start/Stop
                if (Mapping == '5') & (iX == "Start"):
                    df = dfTrimmiX5(df, Span, iX, inside_gene=39, outside_gene=21)
                elif (Mapping == '5') & (iX == "Stop"):
                    df = dfTrimmiX5(df, Span, iX, inside_gene=48, outside_gene=3)
                # todo: correct region for 3' assignemnt
                elif Mapping == '3':
                    df = dfTrimmiX3(df, Span, iX, inside_gene=33, outside_gene=10)
                else:
                    pass


                xt = list(df.index[::6])+[0]

                for i, readLen in enumerate(readLen_l):
                    colors = colorsCheck(colors, readLen)

                    df[readLen].plot(kind='bar', stacked=True, color=colors[readLen], \
                                     width=0.9, legend=True, alpha=a, ax=axes[i])
                    plt.sca(axes[i])
                    axes[i].set_ylabel(Params["Normalised"])

                sns.despine()  # seaborn_aesthetic
                plt.tight_layout()
                fig.savefig(outfig, dpi=300)
            else:
                print("Missing InFile -> {}".format(infile))


def metagPlotsPdf(SRAList, Names, Params):
    #todo: how it works with 3pr mapping
    # using pandas, matplotlib, seaborn, numpy
    makeDirectory("8-MetagPlot")
    sns.set_style("white")    # seaborn_aesthetic
    sns.set_context("paper")  # seaborn_aesthetic

    # Using laTeX to set Helvetica as default font
    from matplotlib import rc
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    ## for Palatino and other serif fonts use:
    # rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)

    from IPython.display import set_matplotlib_formats
    set_matplotlib_formats('pdf', 'svg')

    Span = int(Params["MetagSpan"])
    Mapping = Params["Mapping"]
    dataNorm = Params["Normalised"]  # Mapping 5 or 3 prime end
    rlrange = Params["ReadLenMiN"] + "-" + Params["ReadLenMaX"]  # readlength range -> filename
    readLen_l = [str(i) for i in range(int(Params["ReadLenMiN"]), int(Params["ReadLenMaX"]) + 1)] + ["sum"]

    # colors for plot
    colors = {'25': 'fuchsia', '26': 'blueviolet', '27': 'darkblue', '28': 'b', '29': 'r',
              '30': 'salmon', '31': 'orange', '32': 'olive', '33': 'g', '34': 'tan',
              '35': 'y', 'sum': 'brown', 'ttl': 'brown'
              }

    for iN in Names:
        for iX in ["Start", "Stop"]:
            infile = "7-MetagTbl/" + iN + "_" + Mapping + "-End" + "_" + rlrange + \
                     "_" + dataNorm + "_" + iX + "_iv_Meta_Sum.txt"
            outfig = "8-MetagPlot/" + iN + "-" + Mapping + "-End" + "-" + rlrange + \
                     "-" + dataNorm + "-" + iX + "-iv.pdf"
            outfig_title = "{}".format(iN.replace('_', '-'))

            legend_location = 'upper right' if iX == 'Stop' else 'upper left'

            if os.path.isfile(infile):  # infile exits

                w = 8  # figure width
                h = 1.2 * len(readLen_l)  # figure height

                fig, axes = plt.subplots(nrows=len(readLen_l), figsize=(w, h))

                fig.suptitle(outfig_title, y=0.9, fontsize=12)
                df = pd.read_csv(infile,  index_col=0, sep='\t')
                df.set_index("rel_Pos", inplace=True)

                # Adjust plot for mapping and Start/Stop

                if (Mapping == '5') & (iX == "Start"):
                    df = dfTrimmiX5(df, Span, iX, inside_gene=39, outside_gene=21)
                elif (Mapping == '5') & (iX == "Stop"):
                    df = dfTrimmiX5(df, Span, iX, inside_gene=48, outside_gene=3)
                else: # Mapping == '3'
                    df = dfTrimmiX3(df, Span, iX)

                for i, readLen in enumerate(readLen_l):
                    a = 0.6
                    colors = colorsCheck(colors, readLen)
                    x = df.index
                    y = list(df.loc[:, readLen])
                    axes[i].bar(x, y, color=colors[readLen], alpha=a)
                    axes[i].legend([readLen], loc=legend_location)

                    # colors for guide lines; adjust for beg and end for 5pr
                    b, e = (df.index.min(), df.index.max())
                    # todo: getting axvline colors can be function
                    for k in list(range(b, e+1, 3)):
                        color = 'gray'
                        if k == -12:
                            color = 'g';    a = 0.5
                        elif k == 0:
                            color = 'r';    a = 0.4
                        elif k < 0:
                            color = 'gray'; a = 0.2
                        else:
                            color = 'gray'; a = 0.2

                        # add line after each 3 nt
                        axes[i].axvline(x=k, linewidth=1, alpha=a, color=color)

                    axes[i].set_ylabel(Params["Normalised"])

                sns.despine()  # seaborn_aesthetic
                fig.savefig(outfig, format='pdf', dpi=300, bbox_inches='tight')
                print("{}".format(outfig))
            else:
                print("Missing InFile -> {}".format(infile))
    print("\n")

def read_FASTA(filename, SplitHeader=True):
    """ Reads FastA file and returns a list of tuples, where first
    part is a list of header elements and second seq as a string

    read_FASTA('seqfile.fa', SplitHeader=True)
        [(['gi', '1114567', 'gb', 'NC_00245'],
        'ATATAGGCGCTTGGTGCGCGGCGGGCGCGGCTAGCAGCACCTTTAGTAGCTTTCATCAT'),
        (['gi', '2224567', 'gb', 'NC_22245'],
        'ATTTTTGGGGGGCGCGGCTAGCAGCACCTTTAGTAGCTTTCAAAAAATTTTCAT')]

    info is:
        >gi|1114567|gb|NC_00245

    "Bioinformatics Programming Using Python by Mitchell L Model"
    """
    with open(filename) as file:

        if SplitHeader:
            return [(part[0].split('|'),
                     part[2].replace('\n', ''))
                    for part in
                    [entry.partition('\n')
                     for entry in file.read().split('>')[1:]]]
        else:
            return [(part[0],
                     part[2].replace('\n', ''))
                    for part in
                    [entry.partition('\n')
                     for entry in file.read().split('>')[1:]]]


def read_FASTA_dictionary(filename, SplitHeader=False):
    """ Creates dictionary from FastA file, where key is gi-number and value is seq

    make_indexed_sequence_dictionary('seqfile.fa')
        {'1114567': 'ATATAGGCGCTTGGTGCGCGGCGGGCGCGGCTAGCAGCACCTTTAGTAGCTTTCATCAT',
         '2224567': 'ATTTTTGGGGGGCGCGGCTAGCAGCACCTTTAGTAGCTTTCAAAAAATTTTCAT' }

    read_FASTA by default splits header '|' assuming NCBI entry but
    here read_FASTA do not split header  (SplitHeader=False), i.e. key is the whole name.

    "Bioinformatics Programming Using Python by Mitchell L Model"
    """
    return {info: seq for info, seq in read_FASTA(filename, SplitHeader=False)}


def corrAssignment(SRAList, Names, Params):
    '''
    It reads in interval data _iv adds missing positions, corrects and writes to *.h5
    corrects for offset and writes to new *_assigned_rpm.h5  file
    in this part I keep lines with 0 in table.
    It is not reasonable for human data or bigger genomes
    '''

    makeDirectory("9-Assigncorr")
    makeDirectory("9-Assigncorr/Reports")

    # save_csv= False  #save output to tab delim csv. *.h5 is saved anyway ##?
    rlmin   = int(Params["ReadLenMiN"])
    rlmax   = int(Params["ReadLenMaX"])
    rlrange = str(rlmin) + "-" + str(rlmax)  # read length range 4 filename

    # todo: 3' mapping - not tested - Mapping goes to filename
    Mapping = Params["Mapping"]  # Mapping 5 or 3 prime end
    if Mapping == "3":
        print("\nWARNINGS!\n   3' Mapping is not implemented yet !\n")
        exit(1)
    # todo: Names and column names of OffsetFile must match. Put some check here !
    df = pd.read_csv(Params["OffsetFile"], index_col=0, sep="\t")


    for iN in Names:

        fn_body = iN + "_" + Mapping + "-End_"

        # todo: check *.h5 index status  - how? file exists OR keys structure
        # todo:    a)filename contains "_idx_" b) keys structure "For_raw/V"
        # todo: if not pass it through restructurate_hd5
        # todo: *.h5 2 level index - keys  "/For_rpm/I" ..  def  restructurate_hd5(in.d5, out.h5, close_outfile=True)
        #
        #infile_h5 = "6-AssignRaw/" + fn_body + rlrange + "_iv" + ".h5"
        infile_idx_h5 = "6-AssignRaw/" + fn_body + rlrange + "_idx_iv" + ".h5"
        infile_h5 = infile_idx_h5
        storage = pd.HDFStore(infile_h5, "r")

        readlen_and_offsets = {i: int(df.loc[i, iN]) for i in df[iN].dropna().index}
        rl_l = list(readlen_and_offsets.keys())
        rl_l.sort()  # readlength with periodicity from the table

        fn_body = fn_body + str(min(rl_l)) + "-" + str(max(rl_l))
        #outfile_hdf = "9-Assigncorr/" + fn_body + "_idx_iv" + "_assign_rpm.h5"
        outfile_hdf = "9-Assigncorr/" + fn_body + "_idx_assign_rpm.h5"
        outp_h5 = pd.HDFStore(outfile_hdf, complevel=5, complib="zlib", mode="w")

        LogFileName = "9-Assigncorr/Reports/" + fn_body + "_assign_corr_log.txt"
        LOG_FILE = open(LogFileName, "wt")

        # 2. for Forw & Rev strand searatedly AND for each chr separatedly
        #todo: Normalisation Raw - _rpm_ is hard coded  priority_low
        keys_list = [i for i in storage.keys() if "_rpm" in i]
        keys_for = [i for i in keys_list if "For_" in i]
        keys_rev = [i for i in keys_list if "Rev_" in i]

        # Process Log
        report = "\nInput 1: {}\nRead length included:   {}".format(Params["OffsetFile"], rl_l)
        print(report)
        LOG_FILE.write(report + "\n")

        report = "\nInput 2: {}\nrlmin:   {}\nrlmax:   {}\nName:    {}\nMapping: {}".format(infile_idx_h5, rlmin, rlmax, iN,
                                                                                         Mapping)
        print(report, "\n")
        LOG_FILE.write(report + "\n")

        # 3. get chr length
        # todo: filename for Genome.fa  as parameter
        # todo: 2G problem. Consider use some another FastA format reader. This one opens file and reads it's content
        # todo: Python 3.5 in OSX can't open files bigger thant 2G - problems with human genomes
        #
        genome = read_FASTA_dictionary("0-References/Genome.fa")
        chr_length = {key: len(genome[key]) for key in genome.keys()}

        # columns to include
        columns = ["Chr", "Strand"] + [str(i) for i in rl_l] + ["sum"]  # colum name #'s are str
        read_length_to_use = [str(i) for i in rl_l]

        # 5.Forward str each keys_for
        for key in keys_for:

            Chr = key.split("/")[-1]

            # 5.1 read uncorrected data from h5 to df
            df1 = storage[key]

            # 5.2 reindex
            new_index = list(range(chr_length[Chr]))
            df1 = df1[columns].reindex(new_index)

            # 5.3 apply offset correction
            for rlen in [str(i) for i in rl_l]:
                df1[rlen] = df1[rlen].shift(readlen_and_offsets[int(rlen)])

            # 5.4
            df1["Chr"] = Chr
            df1["Srand"] = "+"
            # works fine for Yeast  :: todo: for bigger genomes like human -> use interval - skip 0 lines
            df1 = df1[read_length_to_use].fillna(0)
            df1.loc[:, 'sum'] = df1.loc[:, read_length_to_use].sum(axis=1)

            # 5.5  write output to h5
            outp_h5[key] = df1

        # Process Log
        report = "Forward done!" ; print(report, "\n"); LOG_FILE.write("\n" + report + "\n")

        # 6. Reverse str each keys_rev
        for key in keys_rev:
            Chr = key.split("/")[-1]

            # 6.1 read uncorrected data from h5 to df
            df1 = storage[key]

            # 6.2 reindex
            new_index = list(range(chr_length[Chr]))
            df1 = df1[columns].reindex(new_index)

            # 6.3 apply offset correction
            for rlen in [str(i) for i in rl_l]:
                df1[rlen] = df1[rlen].shift(-readlen_and_offsets[int(rlen)])

            # 6.4
            df1["Chr"] = Chr
            df1["Srand"] = "+"
            # works fine for Yeast  ::todo human ->iv
            df1 = df1[read_length_to_use].fillna(0)
            df1.loc[:, 'sum'] = df1.loc[:, read_length_to_use].sum(axis=1)

            # 6.5  write output to h5
            outp_h5[key] = df1

        # Process Log
        report = "Reverse done!"; print(report, "\n"); LOG_FILE.write(report + "\n")
        offsets= [readlen_and_offsets[i] for i in rl_l]
        report = "\nOutput: {}\nRead length included:  {}\nOffsets applied:       {}\nName:    {}\nMapping: {}".format(outfile_hdf, \
                                                                                rl_l, offsets, iN, Mapping)

        outp_h5.close()
        storage.close()

        print(report, "\n")
        LOG_FILE.write(report + "\n")
        LOG_FILE.close()


def metagTables2(SRAList, Names, Params):
    """From corrected assingments
    data in h5 file can be interval, i. e. lines with 0 counts are missing
    or continuous -  both index types are possible
    """
    makeDirectory("10-corrMetagTbl")
    makeDirectory("10-corrMetagTbl/Reports")
    # todo: include UTR length req for metagenomic plots priority high
    #utr5 = 30  # length of 5' UTR.  Start codons with UTRs longer than that included
    #utr3 = 30  # length of 3' UTR.  Stop codons with UTRs longer than that included
    # could utr3/5 will replace span later?
    Span = int(Params["MetagSpancorr"])  # nr of nt before and after 5' position of start/stop codons
    #time.sleep(0.1)
    offsettbl = Params["OffsetFile"]

    # todo:assumes 5' :: todo is 3'
    Mapping = "5"    #Params["Mapping"]  # 5 & 3
    # todo: check d I use it default is rpm
    dataNorm = "rpm" #Params["Normalised"]  # "raw" or "rpm"

    LogFileName = "10-corrMetagTbl/Reports/" + "MetagTabl_" + Mapping + "-End_" + "corr_iv_log.txt"
    LOG_FILE = open(LogFileName, "wt")
    #time.sleep(0.1)

    for iN in Names:
        cf1 = cr1 = cf2 = cr2 = 0  # counters
        report = "\nName: {}".format(iN)
        LOG_FILE.write(report + "\n"); print(report)

        # todo: input file name construction relies on offset table. Is it good or bad? What could bemore robust way?
        df = pd.read_csv(offsettbl, index_col=0, sep="\t")
        rl_l    = [i for i in df[iN].dropna().index]  # readlength with periodicity from the table
        fn_body = iN + "_" + Mapping + "-End_"           # filename body
        fn_body+= str(min(rl_l)) + "-" + str(max(rl_l))  #

        # file names
        outf_start = "10-corrMetagTbl/" + fn_body + "_" + dataNorm + "_Start" + "_iv_Meta_Sum.txt"
        outf_stop = "10-corrMetagTbl/" + fn_body + "_" + dataNorm + "_Stop" + "_iv_Meta_Sum.txt"
        infile_h5 = "9-Assigncorr/" + fn_body + "_idx_assign_rpm.h5"
        # corrected assignment
        hd5 = pd.HDFStore(infile_h5, "r")

        # read from table and add 'sum'
        columns = [str(i) for i in df[iN].dropna().index] + ['sum']

        # Empty DataFrames
        meta_start_dff = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)
        meta_start_dfr = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)
        meta_stop_dff  = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)
        meta_stop_dfr  = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)

        report = "Collecting data around genes Start & Stop - Span: {}".format(Span); print(report)
        # Annotation
        tabixfile = pysam.TabixFile("0-References/genome.gtf.gz", parser=pysam.asGTF())

        # Adjust Threshold if for RPM if Normalization is "rpm"
        Threshold = int(Params["MetagThreshold"])  # has to be here
        if dataNorm == "rpm":                      # assumes BAM file is not deleted
            BamName = "5-Aligned/" + iN + ".bam"   # sorted and indexed BAM
            Threshold = raw_metag_threshold_to_rpm(BamName, Threshold)
        else:
            pass

        for Chr in yeastChr():
            key_f = "For_rpm" + "/" + Chr
            key_r = "Rev_rpm" + "/" + Chr
            # It is much faster to read data chr wise in df than access slices from hd5
            # hd5[key].loc[gtf.start - Span:gtf.start + Span, :] # - is slow
            df_f = hd5[key_f]  # Forward
            df_r = hd5[key_r]  # Reverse

            for gtf in tabixfile.fetch(reference=Chr):

                if (gtf.feature == 'stop_codon') & (gtf.strand == '+'):
                    df = pd.DataFrame() # create empty df
                    # enough coverage?
                    if df_f.loc[gtf.start - Span:gtf.start + Span, "sum"].sum() < Threshold:
                        continue
                    else:
                        df = df_f.loc[gtf.start - Span:gtf.start + Span, :]

                    index = range(gtf.start - Span, gtf.start + Span + 1)
                    # def df_framing() makes it interval safe
                    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
                    meta_stop_dff = meta_stop_dff + df  # sum dataframes
                    cf2 += 1

                # Stop codon in Rev strand
                elif (gtf.feature == 'stop_codon') & (gtf.strand == '-'):
                    df = pd.DataFrame()  # create empty df

                    if df_r.loc[gtf.end - Span - 1:gtf.end + Span - 1, "sum"].sum() < Threshold: #-1 correction for rev strand
                        continue
                    else:
                        df = df_r.loc[gtf.end - Span - 1:gtf.end + Span - 1]

                    index = range(gtf.end - Span - 1, gtf.end + Span)  # -1 correction
                    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
                    meta_stop_dfr = meta_stop_dfr + df  # sum dataframes
                    cr2 += 1

                # Sart codon in Forw strand
                elif (gtf.feature == 'start_codon') & (gtf.strand == '+'):
                    df = pd.DataFrame()  # create empty df

                    if df_f.loc[gtf.start - Span:gtf.start + Span, "sum"].sum() < Threshold:
                        continue
                    else:
                        df = df_f.loc[gtf.start - Span:gtf.start + Span, :]

                    index = range(gtf.start - Span, gtf.start + Span + 1)
                    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
                    meta_start_dff = meta_start_dff + df  # sum dataframes
                    cf1 += 1

                # Start codon in Rev strand
                elif (gtf.feature == 'start_codon') & (gtf.strand == '-'):
                    df = pd.DataFrame()  # create empty df

                    if df_r.loc[gtf.end - Span - 1:gtf.end + Span - 1, "sum"].sum() < Threshold:  # Check for data
                        continue
                    else:
                        df = df_r.loc[gtf.end - Span - 1:gtf.end + Span - 1]  # -1 correction for rev strand

                    index = range(gtf.end - Span - 1, gtf.end + Span)  # -1 correction for rev strand
                    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
                    meta_start_dfr = meta_start_dfr + df  # sum dataframes
                    cr1 += 1

                else:
                    pass

        print("Summing up ...")
        LOG_FILE.write("Summing up ...\n")
        # summing up
        meta_start_sum = meta_start_dff + meta_start_dfr
        meta_stop_sum = meta_stop_dff + meta_stop_dfr
        # saving to file
        report = "Around START: {}\nStart sites included  {}".format(outf_start, cf1+cr1);
        LOG_FILE.write(report + "\n"); print(report)
        # print("Sum of saved table: {}".format(int(meta_start_sum["sum"].sum())))
        meta_start_sum.to_csv(outf_start, sep='\t', header=True, index=True)
        meta_start_sum['rel_Pos'] = list(range(-Span, Span + 1))
        meta_start_sum.to_csv(outf_start, sep='\t', header=True, index=True)

        report = "Around STOP:  {}\nStop sites included  {}".format(outf_stop, cf2+cr2)
        LOG_FILE.write(report + "\n"); print(report)
        meta_stop_sum['rel_Pos'] = list(range(-Span, Span + 1))
        meta_stop_sum.to_csv(outf_stop, sep='\t', header=True, index=True)
        print("Done !")
        hd5.close()

    LOG_FILE.close()



Params, SRAs, Names = parseParams("Param.in")
Options             = {1:downloadData, 2: cutAdapt, 3: qualityFilter, 4: ncRNASubtract, 5: genomeAlign, 6:rawAssignment, 
                       7: metagTables, 8: metagPlotspdf, 9:corrAssignment, 10:metagTables2}
# if LaTex not installed use           8: metagPlots   insted  metagPlotspdf

Start               = time.time()

print("\nParameters defined in Param.in:\n")
print_params(Params)
print("\nNames {}\n".format(Names))

for iOpt in range(int(Params["Start"]), int(Params["Stop"]) + 1):
    Options[iOpt](SRAs, Names, Params)
    
    print("Step {} completed! Time taken thus far: {}".format(iOpt, time.time() - Start))
print("\n")
