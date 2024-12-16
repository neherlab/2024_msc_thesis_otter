import pandas as pd
import numpy as np
import click

@click.command()
@click.argument("input_tsv", type=click.File("r"))
@click.argument("dataset_name")
@click.argument("dst")
def format(input_tsv,dataset_name,dst):
    all_data = pd.read_csv(input_tsv,delimiter="\t")
    metadata = all_data[["consensus_gene_name", "consensus_product","average_length(bp)"]].rename({"consensus_gene_name":"gene","consensus_product":"ann","average_length(bp)":"geneLen"},axis="columns")
    locus = all_data.iloc[:,22:]
    max_count = locus.shape[1]
    locuslist=locus.apply(combine_columns,max_count=max_count,axis=1)

    metadata[["locus","count","dupli","sccg","singleton"]]=np.asarray(locuslist.tolist(),dtype="object")
    metadata["geneId"]=metadata.index+1

    metadata=metadata[["geneId","count","dupli","sccg","singleton","locus"]]

    metadata.to_json(f"{dst}/{dataset_name}-pirate.json",orient='records',indent=1)



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
            if ";" in tag:
                for t in tag.split(";"):
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