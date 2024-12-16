import pandas as pd
import numpy as np
import click

@click.command()
@click.argument("input_csv", type=click.File("r"))
@click.argument("dataset_name")
@click.argument("dst")
def format(input_csv,dataset_name,dst):
    all_data = pd.read_csv(input_csv,dtype=str)
    metadata = all_data[["Gene", "Annotation","Avg group size nuc"]].rename({"Gene":"gene","Annotation":"ann","Avg group size nuc":"geneLen"},axis="columns")
    locus = all_data.iloc[:,14:]
    max_count = locus.shape[1]
    locuslist=locus.apply(combine_columns,max_count=max_count,axis=1)

    metadata[["locus","count","dupli","sccg","singleton"]]=np.asarray(locuslist.tolist(),dtype="object")
    metadata["geneId"]=metadata.index+1

    metadata=metadata[["geneId","count","dupli","sccg","singleton","locus"]]

    metadata.to_json(f"{dst}/{dataset_name}-roary.json",orient='records',indent=1)



def combine_columns(row,max_count):
    locus_list = []
    l = 0
    s = 0
    dupli="no"
    sccg="no"
    singleton="no"
    for tag in row:
        if isinstance(tag,str):
            s+=1
            if "\t" in tag:
                for t in tag.split("\t"):
                    locus_list.append(t)
                    l+=1
                dupli="yes"
            else:
                l+=1
                locus_list.append(tag)
    if s==1:
        singleton="yes"
    if (l==max_count) and (dupli=="no"):
        sccg="yes"
    
    return locus_list, l, dupli, sccg, singleton

format()