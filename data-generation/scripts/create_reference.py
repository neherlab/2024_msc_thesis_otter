import argparse
import os
import pandas as pd
from collections import defaultdict

def is_dupli(gene, dupli):
    if gene in dupli:
        return "yes"
    else:
        return "no"
    
def is_singleton(locus):
    g = set()
    for tag in locus:
        g.add(tag.split("_")[1])
    if len(g)==1:
        return "yes"
    else:
        return "no"

def is_sccg(rows,max_count):
    if (rows["count"]==max_count) and (rows["dupli"]=="no"):
        return "yes"
    else:
        return "no"

if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("n_genomes",type=int)
    parser.add_argument("prefix")
    parser.add_argument("combined_file")
    parser.add_argument("tbl_dir")
    parser.add_argument("out_dir")
    args=parser.parse_args()

    print("read combined file:")
    df = pd.read_csv(args.combined_file,header=None,delimiter="\t",names=["gene","tag"])

    # create lookup table for gene to locus_tag
    locus_dict = defaultdict(list)
    for _,item in df.iterrows():
        locus_dict[item.gene].append(item.tag)

    print("lookup table created")

    dupli = []
    tag_dict = {}
    for file in os.listdir(args.tbl_dir):
        print(f"Working on file: {file}")
        df = pd.read_csv(args.tbl_dir+"/"+file,header=None,delimiter="\t",names=["gene","tag"])
        for row in df.itertuples():
            tag_dict[row.tag]=row.gene
        try:
            ds = pd.concat(g for _, g in df.groupby("gene") if len(g) > 1)["gene"]
            for d in ds:
                dupli.append(d)
        except: pass

    dupli = set(dupli)

    print("Create reference")

    reference = pd.DataFrame(columns=["gene","count","dupli","sccg","singleton","locus"])

    reference["gene"]=locus_dict.keys()
    reference["locus"]=locus_dict.values()
    reference["count"]=reference["locus"].apply(len)
    reference["dupli"]=reference["gene"].apply(is_dupli,dupli=dupli)
    reference["singleton"]=reference["locus"].apply(is_singleton)
    reference["sccg"]=reference[["count","dupli"]].apply(is_sccg,max_count=args.n_genomes,axis=1)
    
    fn = args.out_dir+f"/{args.prefix}-reference.json"
    reference.to_json(fn,orient='records',indent=1)
    print(f"Exported reference to {fn}")