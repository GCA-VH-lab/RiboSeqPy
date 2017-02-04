#!/usr/bin/env python
# For Python 2 comment these 2 lines out
#from __future__ import division
#from __future__ import print_function
#
# run code:
#       python Pipeline.py
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

#  Codes backbone from Radhakrishnan, A., et al. Cell (2016)
#! Changes
#! adapted run on Python 3
#! Mapping part is rewritten -> using pysam
#! Annotation (GTF format) is managed by tabix (comes with pysam)
#! Metagenomic tables and plots were added
#! Genome alignment stored in BAM. SamtoBam(iX) *.sam -> sorted & indexed *.bam
#
# This code requires the following programs to be installed on your computer
#
# 1) sra-tools (https://github.com/ncbi/sra-tools/wiki/Downloads)
# 2) cutadapt  (https://cutadapt.readthedocs.io/en/stable/)   #	via anaconda   % conda install cutadapt
# 3) hisat2    (ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads)   version 2.0.5 or higher
# 4) bowtie2   (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
# 5) samtools  (https://github.com/samtools/samtools/)        #	via anaconda   % conda install samtools
# 6) pysam     (https://github.com/pysam-developers/pysam)    # via anaconda   % conda install pysam
#

def cleanFile(File, Condition):

    # gZip a file and delete the un-gZipped version!

    if Condition == "gzip":
        with open(File, 'rb') as FIn, gzip.open(File + ".gz", 'wb') as FOut:
            shutil.copyfileobj(FIn, FOut)

        sp.Popen(["rm", File])

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
    """ returns df what contains values for all positions in the given range
    df1 is condensed df containing positions with values > 0

    :param df1:     condensed df
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


def dfTrimmiX5(df, Span, iX):
    """Truncates Data Frame to fit in figure 5pr """
    x = 3 * 6
    if (Span == 60) & (iX == "Start"):
        df = df[30 - x:90 + x]
        df.reset_index(inplace=True)
        del df['index']
    elif (Span == 60) & (iX == "Stop"):
        df = df[0:60 + 2 * x]
        df.reset_index(inplace=True)
    else:
        pass

    return df


def colorsCheck(dic, key):
    # import numpy.random as rn
    '''Generates random color for missing key'''
    if key not in dic.keys():
        dic[key] = rn.rand(3, 1)
    return dic


def yeastChr():
    # Ordered yeast Chr list short names from ensembl
    return ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','Mito']


def cutAdapt(SRAList, Names, Params):

    makeDirectory("2-Trimmed")
    makeDirectory("2-Trimmed/Reports")

    for iX in Names:
        
        CutAdapt     = ["cutadapt","--discard-untrimmed", "-a","CTGTAGGCACCATCAAT","-o","2-Trimmed/","1-Raw/"]
        CutAdapt[5] += iX + ".fastq.gz"
        CutAdapt[6] += iX + ".fastq.gz"
        
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

        Length = ""
        cutadapt_report = "2-Trimmed/Reports/" + iX + ".txt"
        #avoid error when 2-Trimmed/Reports/*.txt does not exists - happens when start already trimmed input
        if not os.path.exists(cutadapt_report):
            os.system("gzcat  2-Trimmed/" + iX + ".fastq.gz |wc -l > tmp.t")
            lines_in_fastq = int(open('tmp.t', 'r').read()[:-1])
            Length = int(lines_in_fastq / 4)
            print("{} adapters removed elsewhere! \n".format(iX))
        # review: might need to be changed as leads inconsistency in reprot of short/long removed reads nr.
        # caount of all reads in raw fastq file. 2-Trimmed/*.fastq.gz contains reads with adapter - wo. adapter are removed
        elif os.path.exists(cutadapt_report):
            File = open(cutadapt_report)
            Burn = [File.readline() for Idx in range(0, 8)]
            Length = int(Burn[-1][:-1].split(" ")[-1].replace(",", ""))
            File.close()
        else:
            pass  # no options left

        report = "{}   {:>12,}".format(iX, Length)
        LOG_FILE.write(report + "\n"); print(report)

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
                #pass
            elif Len > read_max_len:
                long_reads  += 1
                #pass
            else:
                for IdxL in range(0,Len):
                    Score   = Score*(1 - PHREDDict[PHRED[IdxL]])

                if (Score > Quality):
                    FileOut.write(Identifier + "\n" + Sequence + "\n" + QIdentifier + "\n" + PHRED + "\n")
                else:
                    low_qual +=1

        report  = " Reads len < {}\t   {:>9,}\n".format(read_min_len, short_reads)
        report += " Reads len > {}\t{:>9,}\n".format(read_max_len, long_reads)
        report += " Quality   < {:2.0f}\t{:>9,}".format(-10 * np.log10(1 - Quality), low_qual)
        report += "\n\t\t\t   {:>10,} reads left".format(Length - (short_reads + low_qual))
        print(report, "\n"); LOG_FILE.write(report + "\n\n")

        File.close()
        FileOut.close()


def ncRNASubtract(SRAList, Names, Params):
    
    makeDirectory("4-Subtracted")
    makeDirectory("4-Subtracted/SAM")
    makeDirectory("4-Subtracted/Reports")

    for iX in Names:
        ncRNA    = "0-References/Indexes/ncRNA"
        Input    = "3-Filtered/" + iX + ".fastq"
        Output   = "4-Subtracted/SAM/" + iX + ".sam"
        Unmapped = "4-Subtracted/" + iX + ".fastq"
        
        Bowtie2  = ["bowtie2", "--no-unal", "-p 6", "--un", Unmapped, "-x", ncRNA, "-U", Input,"-S", Output] 
        Subtract = sp.Popen(Bowtie2, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
        Subtract.wait()
        
        FileOut  = open("4-Subtracted/Reports/" + iX + ".txt","w")
        FileOut.write(Subtract.communicate()[1])
        FileOut.close()
        
        cleanFile("3-Filtered/" + iX + ".fastq", Params["Clean"])
        cleanFile("4-Subtracted/SAM/" + iX + ".sam", Params["Clean"])


def genomeAlign(SRAList, Names, Params):

    makeDirectory("5-Aligned")
    makeDirectory("5-Aligned/Reports")

    for iX in Names:
        Genome  = "0-References/Indexes/Genome"
        Input   = "4-Subtracted/" + iX + ".fastq"
        Output  = "5-Aligned/" + iX + ".sam"
        
        Hisat2  = ["hisat2", "--no-unal", "-p 6", "--no-softclip", "--dta", "-x", Genome, "-U", Input, "-S", Output]
        
        Align   = sp.Popen(Hisat2, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)

        FileOut = open("5-Aligned/Reports/" + iX + ".txt", "w")
        FileOut.write(Align.communicate()[1])
        FileOut.close()

        cleanFile("4-Subtracted/" + iX + ".fastq", Params["Clean"])
        
        SamToBam(iX)


def rawAssignment(SRAList, Names, Params):
    # when chromosome names are numbers 1, 2, 3,... it gives warnings when saving PyTable to HDF
    # uncomment when don't want to see warnings  or run    python -W ignore Pipeline_iv.py
    import warnings
    warnings.filterwarnings("ignore")

    makeDirectory("6-AssignRaw")
    makeDirectory("6-AssignRaw/Reports")
    # include_mapped_twice = True/False influence normalization_factor (via that RPM) & mapped reads
    # if True normalisation_factor is computed caounted all mapped reads :: experimental use only
    # if False reads mapped once are counted only :: more conservative and use
    include_mapped_twice = False  # includes reads mapped twice NH:i:2  :: suggested value False
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
        LogFileName = "6-AssignRaw/Reports/" + iN + "_" + Mapping + "-End_" + rlrange + "_iv_log.txt"
        LOG_FILE    = open(LogFileName, "wt")
        # counters for log
        c = c2_sum = c_once = chr_once = total_no = 0
        # empty dataframe for collecting data
        df_for_sum = pd.DataFrame()
        df_rev_sum = pd.DataFrame()
        # Process Log
        report = "\nBamFile: {}\nrlmin:   {}\nrlmax:   {}\nName:    {}\nMapping: {}".format(BamName, rlmin, rlmax, iN,
                                                                                            Mapping)
        LOG_FILE.write(report + "\n"); print(report, "\n")

        # humanChr() gives an ordered list
        for ref in yeastChr():
            c2 = 0
            defF = defaultdict(list)  # DefaultDict  For
            defR = defaultdict(list)  # DefaultDict  Rev
            ForDict = {}  # Collecting  data For
            RevDict = {}  # Collecting  data Rev

            for read in bamfile.fetch(ref):
                chr_once += 1
                if read.get_tag("NH") == 1 :  # NH tag (NH:i:1) tells how many times read are mapped to genome
                    c += 1
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
                elif (read.get_tag("NH") == 2) & include_mapped_twice:
                    c2 += 1
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
                else:
                    pass

            dummy = [0]
            for rlen in range(rlmin, rlmax + 1):
                ForDict[rlen] = Counter(defF.get(rlen, dummy))  # .get() method if rlen
                RevDict[rlen] = Counter(defR.get(rlen, dummy))  # if don't exist use dummy

            df_for = update_df(pd.DataFrame(ForDict), Chr=ref, strand="+")
            df_rev = update_df(pd.DataFrame(RevDict), Chr=ref, strand="-")

            df_for_sum = pd.concat([df_for_sum, df_for], ignore_index=True) # collect summary table
            df_rev_sum = pd.concat([df_rev_sum, df_rev], ignore_index=True) # collect summary table

            # Log_File pr Chr
            report = "{:<5s}\t{:>10,d} reads".format(ref, c)
            LOG_FILE.write(report + "\n"); print(report)
            # Reset/collect counter data
            c_once += c
            c2_sum += c2 # mapped twice
            total_no += chr_once
            c = chr_once = 0

        # Per Name !!!
        # convert int column names to str
        df_for_sum.rename(columns={i: str(i) for i in range(rlmin, rlmax + 1)}, inplace=True)  # num col_names to str
        df_rev_sum.rename(columns={i: str(i) for i in range(rlmin, rlmax + 1)}, inplace=True)  # num col_names to str

        ## Log Report summary
        report = "\nTotal No of reads {:>11,} mapped to genome\n".format(total_no)
        report += "Number   of reads {:>11,d} mapped once to genome\n".format(c_once)
        report += "Number   of reads {:>11,d} mapped twice to genome\n".format(c2_sum)
        report += "Number   of reads {:>11,d} mapped more than counted already\n".format(total_no - (c_once +c2_sum))
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
        if include_mapped_twice:
            l = [0 for read in bamfile.fetch() if read.get_tag("NH") >= 1]  # all mapped reads
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

    LOG_FILE.close()
    bamfile.close()


def metagTables(SRAList, Names, Params):
    makeDirectory("7-MetagTbl")
    makeDirectory("7-MetagTbl/Reports")

    rlmin = int(Params["ReadLenMiN"])
    rlmax = int(Params["ReadLenMaX"])
    Span = int(Params["MetagSpan"])  # nr of nt before and after 5' position of start/stop codons
    Mapping = Params["Mapping"]  # 5 & 3
    dataNorm = Params["Normalised"]  # "raw" or "rpm"
    Threshold = int(Params["MetagThreshold"])

    columns = [str(i) for i in range(rlmin, rlmax + 1)] + ['sum']
    rlrange = str(rlmin) + "-" + str(rlmax)  # readlength range -> filename
    LogFileName = "7-MetagTbl/Reports/" + "MetagTabl_" + Mapping + "-End_" + rlrange + "_iv_log.txt"
    LOG_FILE = open(LogFileName, "wt")

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

        report = "Collecting data around genes Start & Stop - Span: {}".format(Span)
        LOG_FILE.write(report + "\n"); print(report)
        # Annotation
        tabixfile = pysam.TabixFile("0-References/genome.gtf.gz", parser=pysam.asGTF())
        # Adjust Threshold if for RPM if Normalization is "rpm"
        if dataNorm == "rpm":
            BamName = "5-Aligned/" + iN + ".bam"  # sorted and indexed BAM
            Threshold = raw_metag_threshold_to_rpm(BamName, Threshold)
        else:
            pass

        report = "\nMetagThreshold: {:.2f}\t  data normalization {}".format(Threshold, dataNorm)
        LOG_FILE.write(report + "\n"); print(report)

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
        report = "Around START: {}".format(outf_start)
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
    y = 1.8 * len(readLen_l)  # figure height

    for iN in Names:
        for iX in ["Start", "Stop"]:
            infile = "7-MetagTbl/" + iN + "_" + Mapping + "-End" + "_" + rlrange + \
                     "_" + dataNorm + "_" + iX + "_iv_Meta_Sum.txt"
            outfig = "8-MetagPlot/" + iN + "_" + Mapping + "-End" + "_" + rlrange + \
                     "_" + dataNorm + "_" + iX + "_iv.png"
            if os.path.isfile(infile):  # infile exits
                print("\n{}".format(outfig))

                fig, axes = plt.subplots(nrows=len(readLen_l), figsize=(x, y))
                fig.suptitle(outfig, fontsize=14)
                df = pd.DataFrame.from_csv(infile, sep='\t')
                # Adjust plot for mapping and Start/Stop
                if Mapping == '5':
                    df = dfTrimmiX5(df, Span, iX)
                elif Mapping == '3':
                    df = dfTrimmiX3(df, Span, iX)
                else:
                    pass

                ticks = list(df.index[::3])
                labels = list(df.rel_Pos[ticks])

                for i, readLen in enumerate(readLen_l):
                    colors = colorsCheck(colors, readLen)

                    df[readLen].plot(kind='bar', stacked=True, color=colors[readLen], \
                                     width=0.9, legend=True, alpha=a, ax=axes[i])
                    plt.sca(axes[i])
                    plt.xticks(ticks, labels, rotation=0)
                    sns.despine()  # seaborn_aesthetic

                fig.savefig(outfig, dpi=300)
            else:
                print("Missing InFile -> {}".format(infile))


Params, SRAs, Names = parseParams("Param.in")
Options             = {2: cutAdapt, 3: qualityFilter, 4: ncRNASubtract, 5: genomeAlign, 6:rawAssignment, 7: metagTables,
                       8: metagPlots}
Start               = time.time()

print("Params  {}".format(Params))
print("Names {}\n".format(Names))

for iOpt in range(int(Params["Start"]), int(Params["Stop"]) + 1):
    Options[iOpt](SRAs, Names, Params)
    
    print("Step {} completed! Time taken thus far: {}".format(iOpt, time.time() - Start))
print("\n")
