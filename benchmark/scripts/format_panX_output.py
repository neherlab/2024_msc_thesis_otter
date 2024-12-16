import pandas as pd
import click
import json

@click.command()
@click.argument("input_json", type=click.File("r"))
@click.argument("dataset_name")
@click.argument("id_lookup", type=click.File("r"))
@click.argument("dst")
def format(input_json,dataset_name,id_lookup,dst):
    all_data = pd.read_json(input_json)
    max_count = all_data["count"].max()
    out_df = all_data[["geneId","dupli","locus","allGName"]]
    lookup = json.load(id_lookup)

    out_df = out_df.assign(new_locus=out_df["locus"].apply(extract_loci,lookup=lookup))
    out_df.drop(columns=["locus"],inplace=True)
    out_df.rename(columns={"new_locus":"locus"},inplace=True)

    out_df = out_df.assign(count=out_df["locus"].apply(len))
    out_df = out_df.assign(singleton=out_df["locus"].apply(is_singleton))
    out_df = out_df.assign(sccg=out_df[["count","dupli"]].apply(is_sccg,max_count=max_count,axis=1))

    out_df = out_df[["geneId","count","dupli","sccg","singleton","locus"]]
    
    out_df.to_json(f"{dst}/{dataset_name}-panX.json",orient="records",indent=1)


def extract_loci(tags,lookup):
    out = []
    first_tag = next(iter(lookup))
    for tag in tags.split(" "):
        s=None
        c=tag.count("_")
        if c == 1 and lookup.get(tag):
            #case "ta_g"
            s = tag
        elif c == 1:
            #case "strain"_"tag"
            s = tag.split("_")[1]
        elif c == 3:
            #case "strai_n"_"ta_g"
            s = tag.split("_",2)[2]
        elif c == 2 and lookup.get(tag.rsplit("_",1)[1]):
            #case "strai_n"_"tag"
            s = tag.rsplit("_",1)[1]
        elif c == 2 and lookup.get(tag.split("_",1)[1]):
            #case "strain"_"ta_g"
            s=tag.split("_",1)[1]
        elif c==2:
            s = tag
        else:
            #c>2, check for first tag in lookup and split depending on the number of "_"
            d = c-first_tag.count("_")
            s = tag.split("_",d)[d:][0]

        #CHECK IF s IS IN LOOKUP FIRST
        if s and lookup.get(s):
            out.append(lookup.get(s))
        elif s:
            print(f"{s} of {tag} should be in idlookup json")
        else:
            print(f"Can't process {tag}")
    return out

def is_sccg(row,max_count):
    if (row.iloc[0]==max_count) and (row.iloc[1]=="no"):
        return "yes"
    else:
        return "no"
    
def is_singleton(locus):
    g = set()
    for tag in locus:
        g.add(tag.split(".")[0])
    if len(g)==1:
        return "yes"
    else:
        return "no"
    
format()