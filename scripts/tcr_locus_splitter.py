"""
Script to split the GKH TCR library into constant region files.

You'll first need to place the raw fastq files in the raw/ directory and the merged fastq file in the merged/ directory.

Dependencies can be installed using the IgDiscover v1.0.4 conda environment.
"""

import pandas as pd 
import Bio
from Bio import SeqIO
from Bio import Seq
import glob
from tqdm import tqdm
import gzip
import os

base_dir = "../../data/libraries/GKH_TCR/"
raw_dir = "../../data/libraries/GKH_TCR/raw/"
merged_dir = "../../data/libraries/GKH_TCR/merged/"


library_info = pd.DataFrame()
library_info["R1_file"] = glob.glob(os.path.join(raw_dir, "*R1.fastq.gz"))
library_info["loci"] = library_info["R1_file"].apply(lambda x: x.split("_")[-3])
library_info['merged_file'] = library_info["R1_file"].apply(lambda R1: os.path.join(merged_dir, os.path.basename(R1).rsplit("_", 1)[0] + "_merged.fastq.gz"))
library_info['donor'] = library_info['R1_file'].apply(lambda R1: os.path.basename(R1).split("_")[0])

constants = dict()
constants['A'] = "AACCCTGACCCTG"
constants['B'] = "AGGACCTGAAAAACGTG"
constants['D'] = "CAGCCTCATACCA"
constants['G'] = "ACTTGATGCAGATGTT"

def fast2df(fastapath, extra_cols = True):
    """
    Takes a path to a fasta/fastq, returns a pandas dataframe.
    Skips fastq quality indicators
    """
    if "fasta" in fastapath:
        fasttype = "fasta"
    elif "fastq" in fastapath:
        fasttype = "fastq"
    if fastapath[-2:] == "gz":
        return pd.DataFrame([(str(r.description), str(r.seq)) for r in SeqIO.parse(gzip.open(fastapath, "rt"), fasttype)], columns = ['sequence_id', "sequence"])
    else:
        return pd.DataFrame([(str(r.description), str(r.seq)) for r in SeqIO.parse(fastapath, fasttype)], columns = ['sequence_id', "sequence"])

def df2fasta(seq_df, fastapath):
    """
    Takes a dataframe with columns sequence and sequence_id and writes them to the given fasta file.
    """
    records = [Bio.SeqRecord.SeqRecord(seq=Seq.Seq(row["sequence"]), name="", id = row["sequence_id"], description = "") for row in seq_df.to_dict("records")]
    if "fasta" in fastapath:
        fasttype = "fasta"
    elif "fastq" in fastapath:
        fasttype = "fastq"
    if fastapath[-2:] == "gz":
        fastq = gzip.open(fastapath, "wb")
        for record in records:
            fastq.write(str.encode(record.format("fasta")))
        fastq.close()
    else:
        SeqIO.write(records, fastapath, format=fasttype)

def which_constant(seq):
    for constant in constants.keys():
        pos = seq.find(constants[constant])
        if pos != -1:
            return [constant, pos]
    return [None, None]


library_info[["A", "B", "D", "G", "NA"]] = 0
for row in tqdm(library_info.reset_index().to_dict("records")):
    print(row['merged_file'])
    df = fast2df(row['merged_file'])
    # get constant region name and position
    res = df['sequence'].apply(which_constant)
    df['constant_found'] = [el[0] for el in res]
    df['constant_position'] = [el[1] for el in res]
    df['sequence_len'] = df['sequence'].str.len()
    
    # summarize how many of each constant region were found
    constant_ctr = df.groupby(by = 'constant_found', as_index = False).count()[['constant_found', "sequence"]]
    constant_ctr.index = constant_ctr['constant_found']
    constant_ctr = constant_ctr.drop(columns = ['constant_found'])
    constant_ctr = constant_ctr.transpose()
    constant_ctr.at["sequence", "NA"] = int(len(df.query("constant_found.isna()")))
    for col in constant_ctr.columns:
        library_info.loc[row['index'], col] = constant_ctr[col].values[0]
    if row['loci'] == "TCRAD":
        d = os.path.join(base_dir, "TRAD", row['donor'])
        if not os.path.isdir(d):
            os.mkdir(d)
        p = os.path.join(d, f"{row['donor']}_AD.fasta.gz")
        print(p)
        df2fasta(df.query("constant_found == 'A' or constant_found == 'D'").reset_index(drop = True), p)
    else:
        for constant, sub_df in df[df.constant_found.isin(list(row['loci']))].groupby(by = "constant_found"): 
            print(constant)
            d = os.path.join(base_dir, "TR" + constant, row['donor'])
            if not os.path.isdir(d):
                os.mkdir(d)
            p = os.path.join(d, f"{row['donor']}_{constant}.fasta.gz")
            print(p)
            df2fasta(sub_df, p)
    sub_df = df.query("constant_found.isna()")
    d = os.path.join(base_dir, "TRN", row['donor'])
    if not os.path.isdir(d):
        os.mkdir(d)
    p = os.path.join(d, f"{row['donor']}_TRN.fasta.gz")
    print(p)
    df2fasta(sub_df, p)

constant_ctr.to_csv("../data/igdiscover22/GKH_TCR/constant_region_counts.tsv", sep = '\t', index = False)