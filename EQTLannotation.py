#Python 2 
import pandas as pd 
import sys
import urllib2
import urllib
import re
import argparse

def main(args):
    #Following files parse user input file
    file_name=args.Input

    #Open input file and read RSIDs to list
    f=open(file_name,"r")
    rsid_list=[]
    for line in f:
        line=line.strip()
        match = re.search(r"^[r][s][0-9]*$" , line)
        if match:
            rsid_list.append(line)
    f.close()
        
    #GWAS study checker 
    g=open("GWAS","r")
    df = pd.DataFrame(rsid_list)
    GWAS_list=[]
    i=0
    for line in g:
            x=line.split("\t")[21]
            if x in rsid_list:
                GWAS_list.append(x)
    g.close()
    to_append=[]
    for element in rsid_list:
        if element in GWAS_list:
            to_append.append("yes")
        else:
            to_append.append("no")
    df["GWASrecord"]=to_append

    #EXSNPrecord
    #website = urllib2.urlopen("http://www.exsnp.org/data/raw_eQTL.txt")
    #ex_snp={}
    #to_append=[]
    #to_append_2=[]
    #for line in website:
     #   x=("rs"+line.split()[0])
      #  if x in rsid_list:
       #     ex_snp[x]=line.split()[4]
    #for element in rsid_list:
    #    if element in ex_snp.keys():
     #       to_append.append("yes")
      #      to_append_2.append(ex_snp[element])
      #  else:
       #     to_append.append("no")
        #    to_append_2.append("n/a")
    #df["EXSNPrecord"]=to_append
    #df["EXSNPgene"]=to_append_2
    
    #SCAN
    scan={}
    to_append=[]
    to_append_2=[]
    f=open("affy6.dat","r")
    for line in f:
        if line.split()[0] in rsid_list:
            scan[line.split()[0]]=line.split()[3]
    for element in rsid_list:
        if element in scan.keys():
            to_append.append("yes")
            to_append_2.append(scan[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    df["Scan"]=to_append
    df["Scangene"]=to_append_2
    f.close()

    #BloodEQTL CIS
    blood_eqtl={}
    to_append=[]
    to_append_2=[]
    f=open("2012-12-21-CisAssociationsProbeLevelFDR0.5.txt","r")
    for line in f:
        if line.split()[1] in rsid_list:
            if line.split()[1] in blood_eqtl.keys():
                blood_eqtl[line.split()[1]]+=(" ,"+line.split()[13])
            else:
                blood_eqtl[line.split()[1]]=line.split()[13]
    for element in rsid_list:
        if element in blood_eqtl.keys():
            to_append.append("yes")
            to_append_2.append(blood_eqtl[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    df["BloodCisEQTLrecord"]=to_append
    df["BloodCisEQTLgene"]=to_append_2
    f.close()

    #BloodEQTL Trans
    blood_eqtl={}
    to_append=[]
    to_append_2=[]
    f=open("2012-12-21-TransEQTLsFDR0.5.txt","r")
    for line in f:
        if line.split()[1] in rsid_list:
            if line.split()[1] in blood_eqtl.keys():
                blood_eqtl[line.split()[1]]+=(" ,"+line.split()[13])
            else:
                blood_eqtl[line.split()[1]]=line.split()[13]
    for element in rsid_list:
        if element in blood_eqtl.keys():
            to_append.append("yes")
            to_append_2.append(blood_eqtl[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    df["BloodTransEQTLrecord"]=to_append
    df["BloodTransEQTLgene"]=to_append_2
    f.close()

    #GTEX_variant_id_conversion
    f=open("GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt","r")
    variant_id={}
    variant_id_list=[]
    for line in f:
        if line.split()[6] in rsid_list:
            variant_id[line.split()[6]]=(line.split()[2])
    for element in rsid_list:
        if element in variant_id.keys():
            variant_id_list.append(variant_id[element])
        else:
            variant_id_list.append("n/a")
    f.close()

    #GTEX by tissue -lung 
    f=open("Lung.signifpairs.txt","r")
    tissue_eqtl={}
    to_append=[]
    to_append_2=[]
    for line in f:
        if line.split()[0] in variant_id_list:
            if line.split()[0] in tissue_eqtl.keys():
                tissue_eqtl[line.split()[0]]+=(" ,"+line.split()[1])
            else:
                tissue_eqtl[line.split()[0]]=line.split()[1]
    for element in variant_id_list:
        if element in tissue_eqtl.keys():
            to_append.append("yes")
            to_append_2.append(tissue_eqtl[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    df["LungGTEXEQTLrecord"]=to_append
    df["LungGTEXEQTLgene"]=to_append_2
    f.close()
    
    #Adipose_subcutaneous
    f=open("/home/conor/GTEx_Analysis_v7_eQTL/Adipose_Subcutaneous.signifpairs.txt","r")
    tissue_eqtl={}
    to_append=[]
    to_append_2=[]
    for line in f:
        if line.split()[0] in variant_id_list:
            if line.split()[0] in tissue_eqtl.keys():
                tissue_eqtl[line.split()[0]]+=(" ,"+line.split()[1])
            else:
                tissue_eqtl[line.split()[0]]=line.split()[1]
    for element in variant_id_list:
        if element in tissue_eqtl.keys():
            to_append.append("yes")
            to_append_2.append(tissue_eqtl[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    df["Adipose_subGTEXEQTLrecord"]=to_append
    df["Adipose_subGTEXEQTLgene"]=to_append_2
    f.close()
    
    #Cells_transformed_fibroblasts
    f=open("/home/conor/GTEx_Analysis_v7_eQTL/Cells_Transformed_fibroblasts.signifpairs.txt","r")
    tissue_eqtl={}
    to_append=[]
    to_append_2=[]
    for line in f:
        if line.split()[0] in variant_id_list:
            if line.split()[0] in tissue_eqtl.keys():
                tissue_eqtl[line.split()[0]]+=(" ,"+line.split()[1])
            else:
                tissue_eqtl[line.split()[0]]=line.split()[1]
    for element in variant_id_list:
        if element in tissue_eqtl.keys():
            to_append.append("yes")
            to_append_2.append(tissue_eqtl[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    df["FibroblastGTEXEQTLrecord"]=to_append
    df["FibroblastGTEXEQTLgene"]=to_append_2
    f.close()
    
    #Skeletal muscle
    f=open("/home/conor/GTEx_Analysis_v7_eQTL/Muscle_Skeletal.signifpairs.txt","r")
    tissue_eqtl={}
    to_append=[]
    to_append_2=[]
    for line in f:
        if line.split()[0] in variant_id_list:
            if line.split()[0] in tissue_eqtl.keys():
                tissue_eqtl[line.split()[0]]+=(" ,"+line.split()[1])
            else:
                tissue_eqtl[line.split()[0]]=line.split()[1]
    for element in variant_id_list:
        if element in tissue_eqtl.keys():
            to_append.append("yes")
            to_append_2.append(tissue_eqtl[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    df["SkeletalMuscleGTEXEQTLrecord"]=to_append
    df["SkeletalMuscleGTEXEQTLgene"]=to_append_2
    f.close()
    #Skin
    f=open("/home/conor/GTEx_Analysis_v7_eQTL/Skin_Not_Sun_Exposed_Suprapubic.signifpairs.txt","r")
    tissue_eqtl={}
    to_append=[]
    to_append_2=[]
    for line in f:
        if line.split()[0] in variant_id_list:
            if line.split()[0] in tissue_eqtl.keys():
                tissue_eqtl[line.split()[0]]+=(" ,"+line.split()[1])
            else:
                tissue_eqtl[line.split()[0]]=line.split()[1]
    for element in variant_id_list:
        if element in tissue_eqtl.keys():
            to_append.append("yes")
            to_append_2.append(tissue_eqtl[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    df["SkinGTEXEQTLrecord"]=to_append
    df["SkinGTEXEQTLgene"]=to_append_2
    f.close()

    #Blood
    f=open("/home/conor/GTEx_Analysis_v7_eQTL/Whole_Blood.signifpairs.txt","r")
    tissue_eqtl={}
    to_append=[]
    to_append_2=[]
    for line in f:
        if line.split()[0] in variant_id_list:
            if line.split()[0] in tissue_eqtl.keys():
                tissue_eqtl[line.split()[0]]+=(" ,"+line.split()[1])
            else:
                tissue_eqtl[line.split()[0]]=line.split()[1]
    for element in variant_id_list:
        if element in tissue_eqtl.keys():
            to_append.append("yes")
            to_append_2.append(tissue_eqtl[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    df["BloodGTEXEQTLrecord"]=to_append
    df["BloodGTEXEQTLgene"]=to_append_2
    f.close()

    #EQTLs results
    df.to_csv("EQTLresults")

    #TFBS
    f=open("/home/conor/snp2tfbs_JASPAR_CORE_2014_vert.txt","r")
    tfbs_dict={}
    to_append=[]
    to_append_2=[]
    for line in f:
        if line.split()[0] in rsid_list:
            if line.split()[0] in tfbs_dict.keys():
                tfbs_dict[line.split()[0]]+=(" ,"+line.split()[6])
            else:
                tfbs_dict[line.split()[0]]=line.split()[6]
    for element in rsid_list:
        if element in tfbs_dict.keys():
            to_append.append("yes")
            to_append_2.append(tfbs_dict[element])
        else:
            to_append.append("no")
            to_append_2.append("n/a")
    cf = pd.DataFrame(rsid_list)
    cf["SNP2TFBS"]=to_append
    cf["SNP2TFBSgene"]=to_append_2
    f.close()
    cf.to_csv("TFBSresults")
    
    
    
if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("-I","--Input",
                    help="Input RSID file",default="Rnum")
    args = parser.parse_args()
    main(args)

