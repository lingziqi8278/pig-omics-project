#Variance decomposition of gene expression

import os
import limix
import numpy as np
from os.path import join
from pandas import DataFrame
import xarray as xr
import pandas as pd
from numpy import dot


(bim, fam, bed) = limix.io.plink.read("snp_individual",verbose=False)

rna=pd.read_table("rna_for_VC.txt",sep="\t")
position=pd.read_table("gene_position.txt",sep="\t",header=None)
rna_arry=rna.values
pheno = xr.DataArray(rna_arry,[("sample",fam["fid"].tolist()),("gene",position.iloc[:,0].tolist())])
pheno.rename("phenotype")
pheno["chrom"] = ("gene",position.iloc[:,1].tolist())
pheno["start"] = ("gene",position.iloc[:,2].tolist())
pheno["end"] = ("gene",position.iloc[:,3].tolist())

ac=pd.read_table("ac_for_VC.txt",sep="\t")
ac_position=pd.read_table("ac_position.txt",sep="\t",header=None)
ac_arry=ac.values
ac_pheno = xr.DataArray(ac_arry,[("sample",fam["fid"].tolist()),("ac",ac_position.iloc[:,1].tolist())])
ac_pheno.rename("ac")
ac_pheno["chrom"] = ("ac",ac_position.iloc[:,2].tolist())
ac_pheno["start"] = ("ac",ac_position.iloc[:,3].tolist())
ac_pheno["end"] = ("ac",ac_position.iloc[:,4].tolist())

Kh = dot(pheno,pheno.T)
Kh1 = Kh/15509 

G_all=bed.compute()
da = xr.DataArray(G_all.astype('float64'),[("snp",bim["snp"].tolist()),("sample",fam["fid"].tolist())])
da.rename("genotype")
da["chrom"] = ("snp",bim["chrom"].tolist())
da["pos"] = ("snp",bim["pos"].tolist())
da1=da.transpose()
G_all = da1
Kg = dot(G_all, G_all.T)

Y = pheno
AC=ac_pheno
window_size = int(1e6)
vardecs = []

output = open(path+str(name)+'.txt', 'w')
for gene in lysine_group[:]:
        print(".. fit gene "+gene)
        y = Y[:, (Y["gene"] == gene)]
        start = y["start"].item() - window_size
        end = y["end"].item() + window_size
        ac = AC[:, ~(AC["end"] <= start) & ~(AC["start"] >= end) & (AC["chrom"] == int(y["chrom"].item()))]
        if(ac.shape[1]==0):
                continue
        KI = dot(ac, ac.T)    #can be H3K27ac or genetic variants
        vardec = limix.vardec.VarDec(y, "normal")
        vardec.append(KI, "ac")      
        vardec.append(Kg, "geno")
        vardec.append(Kh1, "rna")
        vardec.append_iid("noise")
        vardec.fit(verbose=False)
        vardecs.append(vardec)
output.write(str(vardecs))
output.close()








