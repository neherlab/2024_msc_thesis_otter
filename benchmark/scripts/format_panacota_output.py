import pandas as pd
import click
import json

@click.command()
@click.argument("pangenome_summary", type=click.File("r"))
@click.argument("pangenome", type=click.File("r"))
@click.argument("lookup_json")
@click.argument("dataset_name")
@click.argument("dst")
def format(pangenome_summary,pangenome,lookup_json,dataset_name,dst):
    input_list = pd.read_csv(pangenome_summary)
    metadata=pd.DataFrame()
    locus_dict = {}
    max_count=input_list["sum_quali"].max()

    with open(lookup_json,"r") as f:
        prt2id = json.load(f)    

    with open(pangenome.name) as f:
        clusters = f.readlines()

    for cluster in clusters:
        line = cluster.split(" ")
        locus_list=[]
        #get id from panacota tag
        for tag in line[1:]:
            t = tag.rstrip()
            if t in prt2id:
                locus_list.append(prt2id.get(t))
            else:
                locus_list.append(".".join(t.rsplit("_",1)))
        locus_dict[int(line[0])]=locus_list

    metadata["count"]=input_list["sum_quanti"]
    metadata["dupli"]=input_list["nb_multi"].apply(booleanize)
    metadata["locus"]=input_list["num_fam"].apply(get_tag,dict=locus_dict)
    metadata["geneId"]=metadata.index+1
    metadata["sccg"]=metadata[["count","dupli"]].apply(is_sccg,max_count=max_count,axis=1)
    metadata["singleton"]=input_list["sum_quali"].apply(is_singleton)

    metadata=metadata[["geneId","count","dupli","sccg","singleton","locus"]]

    metadata.to_json(f"{dst}/{dataset_name}-panacota.json",orient='records',indent=1)

def is_sccg(row,max_count):
    if (row[0]==max_count) and (row[1]=="no"):
        return "yes"
    else:
        return "no"

def is_singleton(c):
    if c==1:
        return "yes"
    else:
        return "no"

def booleanize(i):
    if i==0:
        return "no"
    else:
        return "yes"
    
def get_tag(key,dict):
    return dict[key]

format()