import pandas as pd
import numpy as np
import click

@click.command()
@click.argument("input_csv", type=click.File("r"))
@click.argument("dataset_name")
@click.argument("dst")
def format(input_csv,dataset_name,dst):
    all_data = pd.read_csv(input_csv)
    metadata = all_data[["Gene", "Annotation","Avg group size nuc"]].rename({"Gene":"gene","Annotation":"ann","Avg group size nuc":"geneLen"},axis="columns")
    locus = all_data.iloc[:,8:]
    max_count= locus.shape[1]
    locuslist=locus.apply(combine_columns,axis=1)

    metadata[["locus","count","singleton"]]=np.asarray(locuslist.to_list(),dtype="object") #numpy fix
    metadata["geneId"]=metadata.index+1
    metadata["dupli"]=all_data["Avg sequences per isolate"].apply(checkduplicates)
    metadata["sccg"]=metadata.apply(lambda x: is_sccg(x['count'], x['dupli'],max_count), axis=1)
    metadata=metadata[["geneId","count","dupli","sccg","singleton","locus"]]

    metadata.to_json(f"{dst}/{dataset_name}-panta.json",orient='records',indent=1)

def combine_columns(row):
    locus_list = []
    l = 0
    singleton = "no"
    for tag in row:
        if isinstance(tag,str):
            l+=1
            locus_list.append(tag.split("-",2)[-1])  
    if l==1:
        singleton="yes"
    return locus_list, l, singleton

def checkduplicates(n):
    if n>1:
        return "yes"
    else: 
        return "no"
    
def is_sccg(count,dupli,max_count):
    if (count==max_count) and (dupli=="no"):
        return "yes"
    else:
        return "no"


format()