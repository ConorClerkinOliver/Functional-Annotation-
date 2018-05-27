#Python 2 
import pandas as pd 
df=pd.read_csv("EQTLresults")
eqtl=[]
#eqtl+=((df['EXSNPgene']).tolist())
eqtl+=((df['Scangene']).tolist())
eqtl+=((df['BloodCisEQTLgene']).tolist())
eqtl+=((df['BloodTransEQTLgene']).tolist())
eqtl+=((df['LungGTEXEQTLgene']).tolist())
eqtl+=((df['Adipose_subGTEXEQTLgene']).tolist())
eqtl+=((df['FibroblastGTEXEQTLgene']).tolist())
eqtl+=((df['SkeletalMuscleGTEXEQTLgene']).tolist())
eqtl+=((df['SkinGTEXEQTLgene']).tolist())
eqtl+=((df['BloodGTEXEQTLgene']).tolist())
dict_eqtl={}
out=[]
for x in eqtl:
    if ","in str(x):
        for y in str(x).split(","):
            if ":" in str(y):
                out.append(str(y).split(":")[0])                
    if ":" in str(x):
        out.append(str(x).split(":")[0])
    else:
        out.append(str(x).strip())
for q in out:
    if q not in dict_eqtl.keys():
        dict_eqtl[q]=1
    else:
        dict_eqtl[q]+=1

#3dSNP
cf=pd.read_csv("3DSNP")
cf=cf.loc[cf.Tissue.apply(lambda cat: 'Lung' in cat or 'Blood' in cat or "ESC" in cat or "Blood vessel" in cat)]
dict_3d={}
for x in cf.Gene.unique():
    g=cf.loc[cf['Gene'] == x]
    y=(len(g.RSID.unique()))
    dict_3d[x]=y

#TFBSresults"
wf=pd.read_csv("TFBSresults")
TFBS=((wf['SNP2TFBSgene']).tolist())
tfout=[]
dict_tfbs={}
for x in TFBS:
    if ","in str(x):
        for y in str(x).split(","):
            tfout.append(y.strip())
    else:
        tfout.append(str(x).strip())
for q in tfout:
    if q not in dict_tfbs.keys():
        dict_tfbs[q]=1
    else:
        dict_tfbs[q]+=1

w=pd.DataFrame({'EQTLreferences':pd.Series(dict_eqtl),'3DSNP':pd.Series(dict_3d),'TFBSreferences':pd.Series(dict_tfbs)})
w.to_csv("ResultsMeta")
my_genes=w.index.values

#GO enrichment within string 
import sys 
import urllib2
import json

string_api_url = "https://string-db.org/api"
output_format = "json"
method = "enrichment"
species = "9606"
my_app  = "www.aweseome_app.org"

## Construct the request

request_url = string_api_url + "/" + output_format + "/" + method + "?"
request_url += "identifiers=" + "%0d".join(my_genes)
request_url += "&" + "species=" + species
request_url += "&" + "caller_identity=" + my_app

## Call STRING

try:    
    response = urllib2.urlopen(request_url)
except urllib2.HTTPError as err:
    error_message = err.read()
    print error_message
    sys.exit()

## Read and parse the results

result = response.read()
f=open("GOenrichment","w")
if result:
    data = json.loads(result)
    for row in data:
        term = row["term"]
        preferred_names = ",".join(row["preferredNames"])
        fdr = row["fdr"]
        description = row["description"]

        if fdr < 0.01:
            f.write("\t".join([term, preferred_names, str(fdr), description]))

#Construct Network
            
string_api_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "network"
species = "9606"
my_app  = "www.awesome_app.org"

## Construct the request

request_url = string_api_url + "/" + output_format + "/" + method + "?"
request_url += "identifiers=%s" % "%0d".join(my_genes)
request_url += "&" + "species=" + species
request_url += "&" + "caller_identity=" + my_app

try:
    response = urllib2.urlopen(request_url)
except urllib2.HTTPError as err:
    error_message = err.read()
    print error_message
    sys.exit()

## Read and parse the results
g=("Network","w")
line = response.readline()
while line:
    l = line.strip().split("\t")
    p1, p2 = l[2], l[3]
    experimental_score = float(l[10])
    if experimental_score != 0:
        g.write("\t".join([p1,p2, "experimentally confirmed (prob. %.3f)" % experimental_score]))

    line = response.readline()
