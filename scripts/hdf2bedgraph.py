#!/usr/bin/env python
#
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Converts Ribo-Seq coverage to bedgraph format.')
parser.add_argument('-iN', type=str, help='Sample name')
parser.add_argument('-i',  type=str, help='*.hd5 file name')
parser.add_argument('-dtype',  type=str, help='Data type for oufile', default='raw')
parser.add_argument('-col',  type=str, help='column for values - default is "sum"', default='sum')
parser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")
args = parser.parse_args()

if args.verbose:
    print("verbosity turned on")
    
print("\n-iN     prefix:          {}\n-i      input table:     {}\n-dtype  str for outfile: {}\n-col    column:          {}\n".format(args.iN, args.i, args.dtype, args.col))

usage = "./hdf2bedgraph.py -iN WTS1 -i 9-Assigncorr/WTS1_5-End_28-33_idx_assign_rpm.h5 -dtype rpm  -col sum"

if (args.i==None)|(args.iN==None):
     sys.exit("\n  usage:\n\t{}\n".format(usage))

iN = args.iN
infile_h5 = args.i
col = args.col
dtype = args.dtype

storage = pd.HDFStore(infile_h5, "r")

# open file for Forward BedGraph
f1 = iN + '_'+ dtype +'_For.bedgraph'
f2 = iN + '_'+ dtype +'_Rev.bedgraph'
outfile_forward = open(f1, 'w+')
outfile_reverse = open(f2, 'w+')

# chr_list=list(set(df_f.Chr))
chr_list = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV',
            'XV','XVI','Mito']
chr_length = { 'I':230218, 'II':813184, 'III':316620, 'IV':1531933, 'IX':439888, 'Mito':85779, 
'V':576874, 'VI':270161, 'VII':1090940, 'VIII':562643, 'X':745751, 'XI':666816, 'XII':1078177,
'XIII':924431, 'XIV':784333, 'XV':1091291, 'XVI':948066 }
# Generate header
header_1 = 'track type=bedGraph name="{} F {}" description="SacCer3 BedGraph format"\n'.format(iN,dtype)
header_2 = 'track type=bedGraph name="{} R {}" description="SacCer3 BedGraph format"\n'.format(iN,dtype)
# Write header to file
outfile_forward.write(header_1)
outfile_reverse.write(header_2)


def series_2_df(s):
    '''return df with column "sum" '''
    s.name='sum'
    return pd.DataFrame(s)


def df_fill_iv(df1, index, columns):
    """ returns df what contains values for all positions in the given range
    df1 is condensed df containing positions with values > 0
    :param df1:     condensed df, i. e.  don't contain rows with 0 in index
    :param index:   range of genome positions
    :param columns: list of read length + 'sum'
    """
    # create df2
    df2 = pd.DataFrame(0, index=index, columns=columns)
    df1 = df1.add(df2, fill_value=0, axis=1)
    return df1[columns]
    
    
for ref in chr_list: # take one Chr
    print(" {} ...".format(ref))
    
    ######################
    # For Forward strand
    key = "For_rpm/" + ref
    df_ref_f = storage[key]
    # convert series to df
    if type(df_ref_f) == pd.core.series.Series:
        df_ref_f = series_2_df(df_ref_f)
    # For collecting data
    rpm_start = 0    # initial_value
    pos_start = 0    # remember start
    len_incr = 0     # append to length  #! might be omitted
    rpm_current = 0  # current_value     #! might be omitted
    
    # interval -> continious
    index = list(range(chr_length[ref]+1))
    columns = ['sum']
    df_ref_f = df_fill_iv(df_ref_f, index, columns)
    
    # going throug each position for current reference (Chr)
    for i in df_ref_f.index:
        # current value (rpm) in a line 
        rpm_current = df_ref_f.loc[i, col]  
        pos = i       # current position

        # Compare with previous 
        if rpm_start == rpm_current: # same value in a row as before
            len_incr += 1
        else:
            # generate output
            out_line = "{}\t{}\t{}\t{}\n".format(ref, pos_start, pos, rpm_start)
            outfile_forward.write(out_line)
            # redefine variables for next step
            rpm_start = rpm_current
            pos_start = pos
            len_incr  = 0

    ######################
    # For Reverse strand
    key = "Rev_rpm/" + ref
    df_ref_r = storage[key]
    # convert series to df
    if type(df_ref_r) == pd.core.series.Series:
        df_ref_r = series_2_df(df_ref_r)
    # For collecting data
    rpm_start = 0    # initial_value
    pos_start = 0    # remember start
    len_incr = 0     # append to length  #! might be omitted
    rpm_current = 0  # current_value     #! might be omitted
    
    # interval -> continious
    index = list(range(chr_length[ref]+1))
    columns = ['sum']
    df_ref_r = df_fill_iv(df_ref_r, index, columns)
    
    # going throug each position for current reference (Chr)
    for i in df_ref_r.index:
        # current value (rpm) in a line 
        rpm_current = df_ref_r.loc[i, col]  
        pos = i       # current position

        # Compare with previous 
        if rpm_start == rpm_current: # same value in a row as before
            len_incr += 1
        else:
            # generate output
            out_line = "{}\t{}\t{}\t{}\n".format(ref, pos_start, pos, rpm_start)
            outfile_reverse.write(out_line)
            # redefine variables for next step
            rpm_start = rpm_current
            pos_start = pos
            len_incr  = 0
    
print("\n\tAll done!\n")
outfile_forward.close()
outfile_reverse.close()
print("Output BedGrapf files are:\n\t{}\n\t{}\n\n".format(f1,f2))