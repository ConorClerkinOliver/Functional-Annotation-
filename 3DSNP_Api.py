import pandas
import json
import re
import requests
final=pandas.DataFrame({"RSID":[],"Gene":[],"Tissue":[]})
f=open("Rnum","r")
rsid_list=[]
for line in f:
    rsid_list.append(line.strip())
for rsid in rsid_list:
    try:
        response = requests.get("http://cbportal.org/3dsnp/api.do?id="+rsid+"&format=json&type=3dgene")
        data = response.json()
        gene=[]
        rsid_list=[]
        tissue=[]
        try:
            data[0]["data_loop_gene"]
            for x in (data[0]["data_loop_gene"]):
                rsid_list.append(rsid)
                gene.append(x["loopGene"])
                tissue.append(x["loopCellTissue"])
            df=pandas.DataFrame({"RSID":rsid_list,"Gene":gene,"Tissue":tissue})
            final = final.append(df)
        except KeyError:
            continue
    except requests.ConnectionError:
            continue
final.to_csv("3DSNP")
    
