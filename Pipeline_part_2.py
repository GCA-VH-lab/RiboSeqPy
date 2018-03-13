# For Python 2 comment these 2 lines out
#from __future__ import division
#from __future__ import print_function
#
# run code:
#       python Pipeline_part_2.py
#
#       processes steps from correcting assignment (9) up to relative_codon_fold_difference master table (11?)
#

import os
import sys
import time
import gzip
import shutil
import pysam
import subprocess as sp
import numpy.random as rn
import numpy as np
import pandas as pd
import matplotlib as mpl

#todo: add metagene plot for corected and summarised data

def cleanFile(File, Condition):
    # gZip a file and delete the un-gZipped version!

    if Condition == "gzip":
        with open(File, 'rb') as FIn, gzip.open(File + ".gz", 'wb') as FOut:
            shutil.copyfileobj(FIn, FOut)

        sp.Popen(["rm", File])

    if Condition == "bgzip":
        # bgzip No of cores is hard coded
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


def parseParams(Path):
    # Open the parameter file and read in parameters and files to operate on

    File = open(Path)
    SRAList = []
    Names = []
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


def yeastChr():
    # Ordered yeast Chr list short names from ensembl
    return ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','Mito']


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

        report = "\nInput 2: {}\nrlmin:   {}\nrlmax:   {}\nName:    {}\nMapping: {}".format(infile_idx_h5, rlmin, rlmax,
                                                                                            iN, Mapping)
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
    data in h5 file with missing lines, i.e. lines without counts are missing
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
        Threshold = int(Params["MetagThreshold"])  # has to be here # FILTER 3
        report = "Metagenomic threshold {}\n".format(Threshold)
        if dataNorm == "rpm":                      # assumes BAM file is not deleted
            BamName = "5-Aligned/" + iN + ".bam"   # sorted and indexed BAM
            rep = "Estimating normalisation factor from bam file {}".format(BamName); print(rep)
            #todo: feed raw_metag_threshold_to_rpm() with normalisation_factors - reading BAM takse time!
            Threshold = raw_metag_threshold_to_rpm(BamName, Threshold)  # Adjusting FILTER 3  to rpm
        else:
            pass
        report+= "Metagenomic threshold  corrected {}".format(Threshold)
        print(report); LOG_FILE.write(report + "\n")

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
                    if df_f.loc[gtf.start - Span:gtf.start + Span, "sum"].sum() < Threshold: #FILTER 3
                        continue
                    else:
                        df = df_f.loc[gtf.start - Span:gtf.start + Span, :]

                    index = range(gtf.start - Span, gtf.start + Span + 1)
                    # def df_framing() fills missing positions in idex with 0 - makes it interval safe
                    df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
                    meta_stop_dff = meta_stop_dff + df  # sum dataframes
                    cf2 += 1

                # Stop codon in Rev strand
                elif (gtf.feature == 'stop_codon') & (gtf.strand == '-'):
                    df = pd.DataFrame()  # create empty df

                    if df_r.loc[gtf.end - Span - 1:gtf.end + Span - 1, "sum"].sum() < Threshold: #-1 rev str corr #FILTER 3
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

                    if df_f.loc[gtf.start - Span:gtf.start + Span, "sum"].sum() < Threshold: #FILTER 3
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

                    if df_r.loc[gtf.end - Span - 1:gtf.end + Span - 1, "sum"].sum() < Threshold:  #FILTER 3
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


def codonTablesA(SRAList, Names, Params):

    makeDirectory("11-codonTables")
    makeDirectory("11-codonTables/Reports")

    # Input genome and annotation
    genome_file = "0-References/ScerR64-1-1.85.fa"
    annotation = "0-References/genome.gtf.gz"
    logfile_name = "11-codonTables/Reports/codonTables_A.log"

    LOGFILE = open(logfile_name, "w")

    # genome sequence
    genome = {info: seq for info, seq in read_FASTA(genome_file, SplitHeader=False)}

    # Get number of exons for each protein_coding gene
    tabixfile = pysam.TabixFile(annotation, parser=pysam.asGTF())
    missing_keys = [] # check

    Mapping = Params["Mapping"]  # Mapping 5 or 3 prime end
    if Mapping == "3":
        message = "\nWARNINGS!\n   3' Mapping is not implemented yet !\n"
        LOGFILE.write(message + "\n")
        sys.exit(message)

    # Generates dic {gene_id: number_of_exons}
    gene_ids_exon_No = {}

    # todo: include  multi-exon genes
    for ref in yeastChr():
        for gtf in tabixfile.fetch(reference=ref, start=0, end=None):
            if gtf.feature == 'CDS':
                if gtf.gene_id not in gene_ids_exon_No:  # is on in d
                    gene_ids_exon_No[gtf.gene_id] = int(gtf.exon_number)  # pysam gives number as str
                # new exon No is bigger
                if (gtf.gene_id in gene_ids_exon_No) & (gene_ids_exon_No[gtf.gene_id] < int(gtf.exon_number)):
                    gene_ids_exon_No[gtf.gene_id] = int(gtf.exon_number)

    # list of sample names
    ListA = [x.strip(' ') for x in  Params["GroupA"].split(sep=";")]

    for iN  in ListA:
        # read again to avoid half-way iterations
        tabixfile = pysam.TabixFile(annotation, parser=pysam.asGTF())

        fn_body = iN + "_" + Mapping + "-End_"

        # read length from offsets file
        rl_l = readlen_list_from_offset(Params["OffsetFile"], iN)
        fn_body = fn_body + str(min(rl_l)) + "-" + str(max(rl_l))

        # read data in
        infile_h5 = "9-Assigncorr/" + fn_body + "_idx_assign_rpm.h5"
        storage = pd.HDFStore(infile_h5, "r")

        # output file
        outfile_name = "11-codonTables/" + fn_body + "_codon_table_A.txt"
        outfile = open(outfile_name, 'w')

        # ... for library depth
        # enables to specify it  - runs faster - default is 0
        norm_factor = 0  # if 0  estimates from BAM file.
        if norm_factor == 0:
            BamFile = "5-Aligned/" + iN + ".bam"
            print("\n\tProcessing ... {} \n".format(BamFile))
            norm_factor = normalisation_factor_from_bam(BamFile, NHi=1)
        else:  # take a given value
            pass

        # Thresholds
        #todo: instead of giving specific values give 1 value for example "stringency" -- is not transparent what it means!
        # ... for gene from Param.in
        # GeneRpmMeanThr = float(Params["GeneRpmMeanThr"]) # FILTER 1
        # or set dyamic
        GeneRpmMeanThr = 0.3 / norm_factor  # FILTER 1
        # ... for codon
        # mean of raw counts for codon  - 5/3 = 1.666
        # CodonRpmMeanThr =  1.666/norm_factor
        CodonRpmMeanThr = 1.6 / norm_factor  # FILTER 2

        report = "\n{}\nNormalisation factor: {}\nGeneRpmMeanThr: {}\nCodonRpmMeanThr: {}".format(iN,norm_factor,
                                                                                    GeneRpmMeanThr,CodonRpmMeanThr)
        LOGFILE.write(report+"\n"); print(report + "\n")

        ###############################
        # Generate codon based tables #
        ###############################

        # Change header if changing output lines for gene or codon
        header_gene = 'Chr\tExon\tNo_of_exons\tStrand\tGene_id_treated\tgene_1_leftmost\tgene_1_rigthmost\tu_rpm_1'
        header_codon = 'Position_leftmost_1\tcodon_raw_1\tcodon_rpm_1\tcodon_relative_rpm_1\tcodon_P_1\tnorm_factor_1'
        header_codon += '\tamplification_factor_1'
        header = header_gene + "\t" + header_codon + "\n"

        outfile.write(header)

        for ref in yeastChr():
            report = "Chr {:5s} ... ".format(ref); print(report); LOGFILE.write(report+"\n")
            key_f = "For_rpm/" + ref
            key_r = "Rev_rpm/" + ref
            data_f = storage[key_f]
            data_r = storage[key_r]

            for gtf in tabixfile.fetch(reference=ref, start=0, end=None):

                if gtf.feature == 'CDS':
                    #todo: introduce  multiexon genes
                    if gene_ids_exon_No[gtf.gene_id] > 1:  # Skip 2 and more exon genes
                        continue

                    left_most, right_most = gtf.start, gtf.end  # STOP excluded
                    #
                    #left_most, right_most = left_most, right_most - 1  # + strand correction for df.loc[i:j] []-slicing

                    # Retrieve dataframe for gene Check correctness on  "-" strand
                    gene_df = data_f.loc[left_most:right_most-1, ] if gtf.strand == "+" else data_r.loc[left_most:right_most-1, ]

                    if gene_df['sum'].mean() < GeneRpmMeanThr:  # FILTER 1
                        continue

                    # --------------- #
                    # Forward strand
                    if gtf.strand == '+':
                        # get normalization factors for regions
                        collection = gene_region_normalisation_factor(gene_df, left_most, right_most,
                                                                       column='sum', method='mean')

                        norm_factors_coll = tuple([x if x > 0 else np.nan for x in collection])  # repl 0 with nan
                        u_rpm = gene_df["sum"].mean()
                        # todo: consider be consistent taking only CDS and not including STOP codon? left_most+3
                        line_for_gene = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(gtf.contig, gtf.exon_number,
                                                                                str(gene_ids_exon_No[gtf.gene_id]),
                                                                                gtf.strand, gtf.gene_id,
                                                                                left_most, right_most + 3, u_rpm)

                        # For each coding codon on gene, i. e. excluding STOP codon
                        for iPl in range(left_most, right_most, 3):

                            # skip odons with low coverage
                            if gene_df.loc[iPl:iPl + 2,'sum'].mean() < CodonRpmMeanThr:  # FILTER 2
                                continue

                            gene_norm_factor = select_gene_region_norm_factor(norm_factors_coll, iPl, left_most,
                                                                              right_most, strand=gtf.strand)
                            # how much rpms were multiplied to get codon_relative_rpm
                            amplification_factor = 1 / (norm_factor * gene_norm_factor)

                            codon_rpm = gene_df.loc[iPl:iPl + 2, 'sum'].sum()
                            codon_raw = int(codon_rpm * norm_factor) # pseudo raw because input is rpm's original raw is gone
                            codon_relative_rpm = codon_rpm / gene_norm_factor

                            codon_P = genome[ref][iPl:iPl + 3]

                            line_for_codon = "{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}".format(int(iPl), codon_raw,
                                                                                     codon_rpm, codon_relative_rpm,
                                                                                     codon_P, gene_norm_factor,
                                                                                     amplification_factor)

                            line_to_add = line_for_gene + "\t" + line_for_codon + "\n"
                            outfile.write(line_to_add)
                    # --------------
                    # Reverse strand
                    if gtf.strand == '-':
                        # get normalization factors for regions
                        collection = gene_region_normalisation_factor(gene_df, left_most, right_most,
                                                                       column='sum', method='mean')

                        norm_factors_coll = tuple([x if x > 0 else np.nan for x in collection])  # repl 0 with nan
                        u_rpm = gene_df["sum"].mean()
                        # todo: consider be consistent taking only CDS and not including STOP codon? left_most-3
                        line_for_gene = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(gtf.contig, gtf.exon_number,
                                                                                str(gene_ids_exon_No[gtf.gene_id]),
                                                                                gtf.strand, gtf.gene_id,
                                                                                left_most - 3, right_most, u_rpm)

                        # For each coding codon on gene, i. e. excluding STOP codon
                        for iPl in range(left_most, right_most, 3):
                            # skip odons with low coverage
                            if gene_df.loc[iPl:iPl + 2,
                               'sum'].mean() < CodonRpmMeanThr:  # FILTER 2
                                continue

                            gene_norm_factor = select_gene_region_norm_factor(norm_factors_coll, iPl, left_most,
                                                                              right_most, strand=gtf.strand)

                            amplification_factor = 1 / (norm_factor * gene_norm_factor)

                            codon_rpm = gene_df.loc[iPl:iPl + 2, 'sum'].sum()
                            codon_raw = int(codon_rpm * norm_factor)
                            codon_relative_rpm = codon_rpm / gene_norm_factor

                            codon_P = genome[ref][iPl:iPl + 3]
                            codon_P = revcompl(codon_P)
                            line_for_codon = "{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}".format(int(iPl), codon_raw,
                                                                                     codon_rpm, codon_relative_rpm,
                                                                                     codon_P,
                                                                                     gene_norm_factor,
                                                                                     amplification_factor)

                            line_to_add = line_for_gene + "\t" + line_for_codon + "\n"
                            outfile.write(line_to_add)

        report = "\n otufile: {}\n\t{} done!\n".format(outfile_name, iN)
        print(report)
        LOGFILE.write(report)
        outfile.close()

    LOGFILE.close()


def codonTablesB(SRAList, Names, Params):

    makeDirectory("11-codonTables")
    makeDirectory("11-codonTables/Reports")

    # Input genome and annotation
    genome_file = "0-References/ScerR64-1-1.85.fa"
    annotation = "0-References/genome.gtf.gz"
    logfile_name = "11-codonTables/Reports/codonTables_B.log"

    LOGFILE = open(logfile_name, "w")

    # genome sequence
    genome = {info: seq for info, seq in read_FASTA(genome_file, SplitHeader=False)}

    # Get number of exons for each protein_coding gene
    tabixfile = pysam.TabixFile(annotation, parser=pysam.asGTF())
    missing_keys = [] # check

    Mapping = Params["Mapping"]  # Mapping 5 or 3 prime end

    # todo: check 3' mapping
    if Mapping == "3":
        message = "\nWARNINGS!\n   3' Mapping is not implemented yet !\n"
        LOGFILE.write(message + "\n")
        sys.exit(message)

    # todo: include  multi-exon genes
    # Generates dic {gene_id: number_of_exons}
    gene_ids_exon_No = {}
    for ref in yeastChr():
        for gtf in tabixfile.fetch(reference=ref, start=0, end=None):
            if gtf.feature == 'CDS':
                if gtf.gene_id not in gene_ids_exon_No:  # is on in d
                    gene_ids_exon_No[gtf.gene_id] = int(gtf.exon_number)  # pysam gives number as str
                # new exon No is bigger
                if (gtf.gene_id in gene_ids_exon_No) & (gene_ids_exon_No[gtf.gene_id] < int(gtf.exon_number)):
                    gene_ids_exon_No[gtf.gene_id] = int(gtf.exon_number)

    # list of sample names
    ListB = [x.strip(' ') for x in  Params["GroupB"].split(sep=";")]

    for iN  in ListB:
        # read again to avoid half-way iterations
        tabixfile = pysam.TabixFile(annotation, parser=pysam.asGTF())
        # nucleotides before P-site
        nt = int(Params["CodonsBeforePSite"]) * 3
        # whole nucleotide seq length -  part of file name
        sl = str(33 + nt)
        # read length from offsets file
        rl_l = readlen_list_from_offset(Params["OffsetFile"], iN)
        # file name body
        fn_body = iN + "_" + Mapping + "-End_"
        fn_body = fn_body + str(min(rl_l)) + "-" + str(max(rl_l))

        outfile_name = "11-codonTables/" + fn_body + "_codon_table_B_" + sl + ".txt"
        outfile = open(outfile_name, 'w')
        # read data in
        infile_h5 = "9-Assigncorr/" + fn_body + "_idx_assign_rpm.h5"
        storage = pd.HDFStore(infile_h5, "r")

        # codon counters - not used jet
        cP = 0  # plus strand
        cM = 0  # minus strand

        # ... for library depth
        # enables to specify it  - runs faster - default is 0
        norm_factor = 0  # if 0  estimates from BAM file.
        if norm_factor == 0:
            BamFile = "5-Aligned/" + iN + ".bam"
            print("\n\tProcessing ... {} \n".format(BamFile))
            norm_factor = normalisation_factor_from_bam(BamFile, NHi=1)
        else:  # take a given value
            pass

        # Thresholds
        # todo: istead of giving specific values give 1 value for example "stringency" -- is not transparent what it means!
        # ... for gene from Param.in
        # GeneRpmMeanThr = float(Params["GeneRpmMeanThr"]) # FILTER 1
        # or set dyamic
        GeneRpmMeanThr = 0.3 / norm_factor  # FILTER 1
        # ... for codon
        # mean of raw counts for codon  - 5/3 = 1.666
        # CodonRpmMeanThr =  1.666/norm_factor
        CodonRpmMeanThr = 1.6 / norm_factor  # FILTER 2

        report = "\n{}\nNormalisation factor: {}\nGeneRpmMeanThr: {}\nCodonRpmMeanThr: {}".format(iN, norm_factor,
                                                                                                  GeneRpmMeanThr,
                                                                                                  CodonRpmMeanThr)
        LOGFILE.write(report + "\n"); print(report + "\n")

        # Change header if changing output lines for gene or codon
        header_gene = 'Chr\tExon\tNo_of_exons\tStrand\tGene_id_WT\ttranscript_name\tgene_2_leftmost\tgene_2_rigthmost\tu_rpm_2'
        header_codon = 'Position_leftmost_2\tcodon_raw_2\tcodon_rpm_2\tcodon_relative_rpm_2\tcodon_E\tcodon_P\tcodon_A\t'
        header_codon = header_codon + 'sequence\treverse_complement\tpeptide\tnorm_factor_2\tamplification_factor_2'
        header = header_gene + "\t" + header_codon + "\n"

        outfile.write(header)

        for ref in yeastChr():

            report = "Chr {:5s} ... ".format(ref); print(report); LOGFILE.write(report + "\n")
            key_f = "For_rpm/" + ref
            key_r = "Rev_rpm/" + ref
            data_f = storage[key_f]
            data_r = storage[key_r]

            for gtf in tabixfile.fetch(reference=ref, start=0, end=None):

                if gtf.feature == 'CDS':
                    if gene_ids_exon_No[gtf.gene_id] > 1:  # Skip genes containing 2 or more exons
                        continue

                    # --------------- #
                    # Forward strand
                    if gtf.strand == '+':
                        # correct for  '+' strand   df.loc[i:j]  i,j both included
                        left_most, right_most = gtf.start, gtf.end-1

                        gene_df = data_f.loc[left_most:right_most, ]  # Retrieve dataframe for a gene including stop

                        # continue if gene have few RPF's (rpm)
                        if gene_df['sum'].mean() < GeneRpmMeanThr:  # FILTER 1
                            continue

                        # get normalization factors for regions
                        collection = gene_region_normalisation_factor(gene_df, left_most, right_most,
                                                                       column='sum', method='mean')
                        # replace 0 with NaNs to avoid div by 0
                        norm_factors_coll = tuple([x if x > 0 else np.nan for x in collection])

                        # --------------- #
                        #  Line_for_gene  #
                        #### Chr; Exon; No_of_exons; Strand; Gene_id_treated; gleft_most; rigth_most
                        # Changing here update codonBased_df columns !!!
                        u_rpm = gene_df["sum"].mean()
                        line_for_gene = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(gtf.contig, gtf.exon_number,
                                                                                    str(gene_ids_exon_No[gtf.gene_id]),
                                                                                    gtf.strand, gtf.gene_id,
                                                                                    gtf.transcript_name,
                                                                                    left_most, right_most+3, u_rpm)
                        # -------------------------------------------------------- #
                        # For each coding codon on gene, i. e. excluding STOP codon
                        for iPl in range(left_most, right_most, 3):

                            # skipp codons with low coverage
                            if gene_df.loc[iPl:iPl + 2,'sum'].mean() < CodonRpmMeanThr:  # FILTER 2
                                continue

                            cP += 1  # counting codons
                            gene_norm_factor = select_gene_region_norm_factor(norm_factors_coll, iPl,
                                                                              left_most, right_most)
                            amplification_factor = 1 / (norm_factor * gene_norm_factor)

                            # calculate
                            codon_rpm = gene_df.loc[iPl:iPl + 2, 'sum'].sum()
                            codon_raw = int(codon_rpm * norm_factor)
                            codon_relative_rpm = codon_rpm/gene_norm_factor

                            # Special case when E-site is empty
                            codon_E = np.nan if iPl == left_most else genome[ref][iPl - 3:iPl]
                            codon_P = genome[ref][iPl:iPl + 3]      # P site
                            codon_A = genome[ref][iPl + 3:iPl + 6]  # A site

                            pos_30_codons_back = iPl - 30

                            if left_most > pos_30_codons_back:
                                pos_30_codons_back = left_most

                            sequence = genome[ref][pos_30_codons_back:iPl + 3 + nt]  # nt nucleotides before P site
                            revcomp = np.nan  # + strand reverse complement is empty
                            translation = translate_DNA(sequence)

                            # ---------------- #
                            #  Line_for_codon  #
                            line_for_codon = "{}\t{}\t{}\t{}\t".format(int(iPl), codon_raw, codon_rpm,
                                                                       codon_relative_rpm)
                            line_for_codon += "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(codon_E, codon_P, codon_A, sequence,
                                                                                  revcomp, translation, gene_norm_factor)
                            line_for_codon += "\t{:.3f}".format(amplification_factor)
                            line_to_add = line_for_gene + "\t" + line_for_codon + "\n"
                            outfile.write(line_to_add)

                    if gtf.strand == '-':
                        # left_most, right_most = left_most, right_most+1
                        # no correction for  '-' strand
                        left_most, right_most = gtf.start, gtf.end

                        gene_df = data_r.loc[left_most:right_most, ]

                        # get normalization factors for regions (beg, body, end)
                        # select_gene_region_norm_factor() is strand aware and flips beg and end to fit "-" strand
                        collection = gene_region_normalisation_factor(gene_df, left_most, right_most,
                                                                       column='sum', method='mean')

                        norm_factors_coll = tuple([x if x > 0 else np.nan for x in collection])  # repl 0 with nan
                        u_rpm = gene_df["sum"].mean()
                        # todo: consider be consistent taking only CDS and not including STOP codon? left_most-3
                        line_for_gene = "{}\t{}\t{}\t{}\t{}\t{}".format(gtf.contig, gtf.exon_number,
                                                                        str(gene_ids_exon_No[gtf.gene_id]),
                                                                        gtf.strand, gtf.gene_id, gtf.transcript_name)

                        line_for_gene += "\t{}\t{}\t{}".format(left_most - 3, right_most, u_rpm)

                        # For each coding codon on gene, i. e. excluding STOP codon
                        for iPl in range(left_most, right_most, 3):
                            # skip odons with low coverage
                            if gene_df.loc[iPl:iPl + 2, 'sum'].mean() < CodonRpmMeanThr:  # FILTER 2
                                continue

                            cM += 1  # counting codons
                            gene_norm_factor = select_gene_region_norm_factor(norm_factors_coll, iPl, left_most,
                                                                              right_most, strand=gtf.strand)

                            amplification_factor = 1 / (norm_factor * gene_norm_factor)

                            codon_rpm = gene_df.loc[iPl:iPl + 2, 'sum'].sum()
                            codon_raw = int(codon_rpm * norm_factor)
                            codon_relative_rpm = codon_rpm / gene_norm_factor

                            line_for_codon = "{}\t{}\t{}\t{}\t".format(int(iPl), codon_raw, codon_rpm, codon_relative_rpm)
                            # Get sequences
                            codon_E = np.nan if iPl + 3 == right_most else revcompl(genome[ref][iPl + 3:iPl + 6])  # E site
                            codon_P = revcompl(genome[ref][iPl:iPl + 3])  # P site
                            codon_A = revcompl(genome[ref][iPl - 3:iPl])  # A site

                            pos_30_codons_forward = iPl + 30

                            pos_30_codons_forward = right_most if pos_30_codons_forward > right_most else iPl + 30
                            #
                            sequence = genome[ref][iPl - nt:pos_30_codons_forward + 3]  # 11 codons including P-site & nt before P-Site
                            seqrevcomp = revcompl(sequence)
                            translation = translate_DNA(seqrevcomp)

                            # ---------------- #
                            #  Line_for_codon  #
                            line_for_codon += "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(codon_E, codon_P, codon_A, sequence,
                                                                                  seqrevcomp, translation, gene_norm_factor)
                            line_for_codon += "\t{:.3f}".format(amplification_factor)

                            line_to_add = line_for_gene + "\t" + line_for_codon + "\n"
                            outfile.write(line_to_add)

        report = "\n otufile: {}\n\t{} done!\n".format(outfile_name, iN)
        print(report)
        LOGFILE.write(report)
        outfile.close()

    LOGFILE.close()


def codonTablesAB(SRAList, Names, Params):
    # running functions for  tables A and B

    codonTablesA(SRAList, Names, Params)
    codonTablesB(SRAList, Names, Params)


def masterTable(SRAList, Names, Params):
    #
    makeDirectory("12-MasterTables")
    makeDirectory("12-MasterTables/Reports")

    Mapping = Params["Mapping"]
    nt = int(Params["CodonsBeforePSite"]) * 3
    sl = str(33 + nt) # whole nucleotide seq length -  part of file name

    # list of sample names
    ListA = [x.strip(' ') for x in Params["GroupA"].split(sep=";")]
    ListB = [x.strip(' ') for x in Params["GroupB"].split(sep=";")]

    for i, iNa in enumerate(ListA):
        iNb   = ListB[i]

        # inut table names
        # A
        fn_bodyA  = iNa + "_" + Mapping + "-End_"
        rl_l      = readlen_list_from_offset(Params["OffsetFile"], iNa) # read length from offsets file
        fn_bodyA += str(min(rl_l)) + "-" + str(max(rl_l))
        # B
        fn_bodyB  = iNb + "_" + Mapping + "-End_"
        rl_l      = readlen_list_from_offset(Params["OffsetFile"], iNb)
        fn_bodyB += str(min(rl_l)) + "-" + str(max(rl_l))
        # input table names
        infile_B = "11-codonTables/" + fn_bodyB + "_codon_table_B"+ "_" + sl + ".txt"
        infile_A = "11-codonTables/" + fn_bodyA + "_codon_table_A" + ".txt"
        report = "Replica No {}\n  TReated: {}\nunTreated: {}".format(i + 1, iNa, iNb)
        print(report)
        # output file
        outfile_csv        = "12-MasterTables/Master_table_" + iNa + "_" + iNb + ".txt"
        outfile_cleaned_csv= "12-MasterTables/Master_table_" + iNa + "_" + iNb + "_cleaned.txt"
        #outfile            = open(outfile_csv, 'w')
        #outfile_cleaned    = open(outfile_cleaned_csv, 'w')
        # List of columns to be included into Master Table
        columns_needed = ['Chr', 'Exon', 'No_of_exons', 'Strand', 'Gene_id_treated', 'codon_relative_fd',
                      'gene_1_leftmost', 'gene_1_rigthmost', 'Position_leftmost_1', 'codon_raw_1',
                      'codon_rpm_1', 'codon_relative_rpm_1', 'u_rpm_1', 'norm_factor_1',
                      'Gene_id_WT', 'transcript_name', 'Position_leftmost_2', 'codon_raw_2', 'codon_rpm_2',
                      'codon_relative_rpm_2', 'u_rpm_2', 'norm_factor_2', 'codon_E', 'codon_P', 'codon_P_1',
                      'codon_A', 'sequence', 'reverse_complement', 'peptide', 'amplification_factor_1',
                      'amplification_factor_2']

        # Columns convert to int
        col2int = ['Exon', 'No_of_exons', 'gene_1_leftmost', 'gene_1_rigthmost', 'Position_leftmost_1',
               'Position_leftmost_2']

        print("Reading input")
        # Read files in. produced in step 5_
        print("unTreated: {}".format(infile_B))
        df_wt = pd.read_table(infile_B)
        print("  Treated: {}".format(infile_A))
        df_tr = pd.read_table(infile_A)
        print("..")
        # replace 0 with NaN's
        columns = ['codon_raw_2', 'codon_rpm_2', 'codon_relative_rpm_2']
        for col in columns:
            df_wt.loc[df_wt[col] == 0, col] = np.nan

        columns = ['codon_raw_1', 'codon_rpm_1', 'codon_relative_rpm_1']
        for col in columns:
            df_tr.loc[df_tr[col] == 0, col] = np.nan

        # reduce dataframes a bit
        df_wt = df_wt.dropna(subset=['codon_relative_rpm_2']).reset_index()
        df_tr = df_tr.dropna(subset=['codon_relative_rpm_1']).reset_index()

        # Later use for index
        df_tr['Position'] = df_tr['Position_leftmost_1']
        df_wt['Position'] = df_wt['Position_leftmost_2']
        # Create DF for collectiong subframes
        final_df = pd.DataFrame(columns=columns_needed)

        for ref in yeastChr():  # CHANGE IT
            print("{:>6s}".format(ref))
            # split to Frow & Rev because in some cases Position No are overlaping
            # Duplicated Positions appear because of overlapping gene annotation in the same frame
            # For example Genes YBL069W & YBL068W-A  in Chr 'II'
            # it's rare and don't influence owerall result

            # R# Reverse strand operationas are commented out
            dfA_Ftr = df_tr[(df_tr.Chr == ref) & (df_tr.Strand == '+')].drop_duplicates(subset='Position',
                                                                                        keep='last').set_index('Position')
            dfA_Rtr = df_tr[(df_tr.Chr == ref) & (df_tr.Strand == '-')].drop_duplicates(subset='Position',
                                                                                        keep='last').set_index('Position')
            dfB_Fwt = df_wt[(df_wt.Chr == ref) & (df_wt.Strand == '+')].drop_duplicates(subset='Position',
                                                                                        keep='last').set_index('Position')
            dfB_Rwt = df_wt[(df_wt.Chr == ref) & (df_wt.Strand == '-')].drop_duplicates(subset='Position',
                                                                                        keep='last').set_index('Position')
            # Some columns are duplicated, some are uniqe in WT. take what is needed
            columns_wt = ['Gene_id_WT', 'transcript_name', 'Position_leftmost_2', 'codon_raw_2', 'codon_rpm_2',
                          'codon_relative_rpm_2', 'u_rpm_2', 'norm_factor_2', 'codon_E', 'codon_P', 'codon_A',
                          'sequence', 'reverse_complement', 'peptide', 'amplification_factor_2']

            frame_F = [dfA_Ftr, dfB_Fwt[columns_wt]]
            frame_R = [dfA_Rtr, dfB_Rwt[columns_wt]]
            result_F = pd.concat(frame_F, axis=1, join='inner')
            result_R = pd.concat(frame_R, axis=1, join='inner')
            result_F['codon_relative_fd'] = result_F.codon_relative_rpm_1 / result_F.codon_relative_rpm_2
            result_R['codon_relative_fd'] = result_R.codon_relative_rpm_1 / result_R.codon_relative_rpm_2

            result_F.reset_index(inplace=True)
            result_R.reset_index(inplace=True)

            final_df = pd.concat([final_df, result_F[columns_needed]], ignore_index=True)
            final_df = pd.concat([final_df, result_R[columns_needed]], ignore_index=True)

        print("\nSome Final adjustments ...")
        # convert some float to int
        final_df[col2int] = final_df[col2int].astype(int)
        # renme some columns
        # df.rename(index=str, columns={"A": "a", "B": "c"})
        columns = {'gene_1_leftmost': 'gene_leftmost', 'gene_1_rigthmost': 'gene_rightmost'}
        final_df.rename(columns=columns, inplace=True)

        final_df.sort_values('codon_relative_fd', ascending=False, inplace=True)
        # final_df.to_excel("6_Master_table_WT_met_S1S2.xlsx")
        print("\nSaving results ...")
        final_df.to_csv(outfile_csv, sep='\t')
        report = "\t{}\nCodons included {:>10,}".format(outfile_cleaned_csv, final_df.shape[0])
        print(report)

        ## do some cleaning
        # using mask to filter out specific len of amino acid sequences
        mask = (final_df.peptide.str.contains('_') == False)
        final_df = final_df.loc[mask]
        # remove shorter than max len
        # 33 + nt - nt determines codons before P-Site  int(Params["CodonsBeforePSite"])*3
        mask = (final_df.sequence.str.len() == 33 + nt)
        final_df = final_df.loc[mask]

        final_df.to_csv(outfile_cleaned_csv, sep='\t')

        report = "\t{}\nCodons included {:>10,}".format(outfile_cleaned_csv, final_df.shape[0])
        print(report)

    print("\nDone!\n")

#---- START helper functions ----#
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


def raw_metag_threshold_to_rpm(BamName, Threshold):
    """
    Converts raw MetagThresold to rpm MetagThreshold

    :param BamName:    BAM file used for caclulating normalization factor
    :param Threshold:  raw threshold
    :return: normalized threshold to fit with rpm normlized data
    todo: check file for normalisation factors. it is less burden especially with big bam files
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


def dfTrimmiX5(df, Span, iX, c=10):
    """Truncates Data Frame to fit in figure 5pr """
    if iX == "Start":
        return  df.loc[-24:Span-c, ]
    elif iX == "Stop":
        return  df.loc[-Span:24-c, ]
    else:
        return df


def dfTrimmiX3(df, Span, iX):
    """Truncates Data Frame to fit in figure 3pr"""
    #todo: rewrite it to accept "rel_Pos" index similarly to dfTrimmiX5()
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


def translate_DNA_codon(codon):
    ''' Translates DNA triplet to amino acid (20) by using single
    lettr code where underline "_" corresponds to STOP codon.
    '''
    return DNA_codon_table[codon]


def aa_generator_DNA(dnaseq):
    """Return a generator object that produces an amino acid by translating
    the next three characters of dnaseq each time next is called on it"""
    return (translate_DNA_codon(dnaseq[n:n+3])
            for n in range(0, len(dnaseq), 3))


def translate_DNA(dnaseq):
    """Translate dnaseq into amino acid symbols"""

    gen = aa_generator_DNA(dnaseq)
    seq = ''
    aa = next(gen, None)
    while aa:
        seq += aa
        aa = next(gen, None)
    return seq


DNA_codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}


def complement(seq):
    ''' Returns complament sequence of DNA
        Expects DNA string as input
        def dictionary is  extended IUPAC
    '''
    basecomplement = {'A':'T','C':'G','G':'C','T':'A','Y':'R','R':'Y','M':'K','K':'M','W':'W','V':'B','B':'V','H':'D','D':'H','N':'N',
                      'a':'t','c':'g','g':'c','t':'a','y':'r','r':'y','m':'k','k':'m','w':'w','v':'b','b':'v','h':'d','d':'h','n':'n' }
    return ''.join([basecomplement[base] for base in seq])


def revcompl(seq):
    ''' returns inverted string of complament DNA
    '''
    return complement(seq[::-1])


def gene_region_normalisation_factor(df, gene_beg, gene_end, column='sum', method='mean'):
    '''
    Returns tuple of noramlisation factors for CDS (begin_region, body_region, end_region).

    This method (df.loc[i:j]) use closed end slicing [], i.e. both ends i & j are included

    Normalisation factor calculated basd on a given method over three regions.

    Valid methods are:
        'mean', 'sum' & 'max'

    Regions are:
        begin_region - 10 codons from begin
        end_region   - 10 codons from end
        body_region  - between begin_ - & end_region

    Work for Forward strand genes.

   For Reverse strand returned tuple have to reverted afterwards.

          norm_f_beg, norm_f_body, norm_f_end = norm_f_end, norm_f_body, norm_f_beg

    It is a part of Riob-Seq data conversion.
    '''

    if method == 'mean':
        return df.loc[gene_beg:gene_beg + 30, column].mean(), df.loc[gene_beg + 30:gene_end - 30, \
                                                              column].mean(), df.loc[gene_end - 30:gene_end,
                                                                              column].mean()

    elif method == 'max':
        return df.loc[gene_beg:gene_beg + 30, column].max(), df.loc[gene_beg + 30:gene_end - 30, \
                                                             column].max(), df.loc[gene_end - 30:gene_end, column].max()

    elif method == 'sum':
        return df.loc[gene_beg:gene_beg + 30, column].sum(), df.loc[gene_beg + 30:gene_end - 30, \
                                                             column].sum(), df.loc[gene_end - 30:gene_end, column].sum()
    else:
        print("method_Error in function gene_region_normalisation_factor()!!!\nNo proter method was provided.")
        print("Valid methods are: 'mean', 'sum' & 'max'\n")


def select_gene_region_norm_factor(norm_factors_coll, Position_leftmost, left_most, right_most, strand="+"):
    """
    norm_factors_coll:  tuple of normalisation factors for gene regions"""

    if len(norm_factors_coll) != 3:
        print("ERROR! Expect 3 fnormalisation factors but give ins {}\t{}".format(len(norm_factors_coll)))
        sys.exit(1)

    if Position_leftmost <= left_most + 30:  # '+' beginning - '-' end
        return norm_factors_coll[0] if strand == "+" else norm_factors_coll[2]
    elif (Position_leftmost > left_most + 30) & (Position_leftmost < right_most - 30):  # main body
        return norm_factors_coll[1]
    elif Position_leftmost >= right_most - 30:  # '+' end - '-' beginning
        return norm_factors_coll[2] if strand == "+" else norm_factors_coll[0]
    else:  # it can't hapen but anyway
        print("ERROR! Position {} is out of range {}\t{}".format(Position_leftmost, left_most, right_most + 3))


def gene_normalisation_factor(df, gene_beg, gene_end, column='sum', method='mean'):
    '''
    Returns tuple of noramlisation factors for CDS.

    This method (df.loc[i:j]) use closed end slicing [], i.e. both ends i & j are included

    Normalisation factor calculated basd on a given method over three regions.

    Valid methods are:
        'mean', 'sum' & 'max'

    It is a part of Riob-Seq data conversion.
    '''

    if method == 'mean':
        return df.loc[gene_beg:gene_end, column].mean()

    elif method == 'median':
        return df.loc[gene_beg:gene_end, column].median()

    elif method == 'sum':
        return df.loc[gene_beg:gene_end, column].sum()
    else:
        print("method_Error in function gene_region_normalisation_factor()!!!\nNo proter method was provided.")
        print("Valid methods are: 'mean', 'sum' & 'median'\n")


def reads_count_in_bam(BamName, NHi=1):
    bamfile = pysam.AlignmentFile(BamName, "rb")  # open BAM file
    c = 0
    for read in bamfile.fetch():
        if read.get_tag("NH") <= NHi:  # mapped once if 1; mapped once and twice if 2
            c += 1
    return c


def normalisation_factor_from_bam(BamName, NHi=1):
    return reads_count_in_bam(BamName, NHi) / (10 ** 6)


def readlen_list_from_offset(OffsetFile, iN):
    '''returns list of readlength from offset files for specified sample'''

    df = pd.read_csv(Params["OffsetFile"], index_col=0, sep="\t")
    readlen_and_offsets = {i: int(df.loc[i, iN]) for i in df[iN].dropna().index}

    return list(readlen_and_offsets.keys())


Params, SRAs, Names = parseParams("Param.in")
Options = {9: corrAssignment, 10: metagTables2, 11: codonTablesAB, 12: masterTable}
Start = time.time()

print("Params  {}".format(Params))
print("Names {}\n".format(Names))

#codonTablesA(SRAList, Names, Params)

for iOpt in range(int(Params["Start"]), int(Params["Stop"]) + 1):

    Options[iOpt](SRAs, Names, Params)

    print("Step {} completed! Time taken thus far: {}".format(iOpt, time.time() - Start))
print("\n")
