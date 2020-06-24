#! /usr/bin/env/ python

import re
import sys
import pandas as pd

vcfFile = sys.argv[1]
avFile = vcfFile + ".avinput"
anFile = vcfFile + ".intervar.hg19_multianno.txt"
inFile = vcfFile + ".intervar.hg19_multianno.txt.intervar"

vcf = pd.read_csv(vcfFile, sep="\t", comment="#", header=None)
av_tmp = pd.read_csv(avFile, sep="\t", header=None)
annovar = pd.read_csv(anFile, sep="\t",skiprows=1, header=None)
intervar_tmp = pd.read_csv(inFile, sep="\t")

with open(anFile, "r") as f:
    head = [next(f) for x in range(2)]

n_missing = len(head[1].strip().split("\t")) - len(head[0].strip().split("\t"))
annovar_cols =  head[0].strip().split("\t") + ["OtherInfo" + str(i) for i in range(n_missing)]

annovar.columns = annovar_cols
intervar = intervar_tmp.iloc[:,[0,1,2,13]]
intervar.columns = ["Chr", "Start", "End", "Intervar"]
del intervar_tmp

avinput = av_tmp.iloc[:,range(3)]
avinput.columns = ["Chr", "Start", "End"]
vcf.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
vcf_ext = pd.concat([avinput, vcf], axis=1)
del av_tmp

def acmg(s):
  pattern = "InterVar: (.*?) PVS1"
  m = re.search(pattern, s).group(1)
  return(m)

intervar["acmg_intervar"] = intervar.Intervar.apply(acmg)
intervar.drop_duplicates(inplace=True)
intervar.drop(labels=["Intervar"], inplace=True, axis="columns")

out = pd.merge(vcf_ext, pd.merge(intervar, annovar))

out.to_csv( vcfFile + "_annotations.tsv", sep="\t", index=False)
