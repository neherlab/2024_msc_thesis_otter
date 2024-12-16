import argparse
import os
import pandas as pd
import numpy as np
import json

import matplotlib
import matplotlib.pyplot as plt


# --------------------------- #
#          methods
# --------------------------- #

#hash a list of locus tags
def get_hash(locus_list):
    locus_list.sort()
    return hash(frozenset(locus_list))

#returns a set without nan values
def set_wo_nan(geneId_list):
    s = set(geneId_list)
    #remove nan values as nan!=nan
    return {x for x in s if x==x}

#make figure title based on filter
def make_title(ref,dataset,f):
    if f:
        # filter by lower and upper bounds (inclusive)
        # if only one value is specified use clusters up to hi or from lo on
        if isinstance(f,tuple):
            lo = f[0]
            hi = f[1]
            if hi and lo:
                return f"distribution of {ref} clusters ({lo}<=size<={hi})\ndataset {dataset}"
            elif hi:
                return f"distribution of {ref} clusters (size<={hi})\ndataset {dataset}"
            elif lo:
                return f"distribution of {ref} clusters ({lo}<=size)\ndataset {dataset}"
        # filter by cluster size
        elif isinstance(f,int):
            return f"distribution of {ref} clusters (size={f})\ndataset {dataset}"
        # filter by sccg or singleton (column specifier)
        elif isinstance(f,str):
           return f"distribution of {ref} clusters ({f})\ndataset {dataset}"
    else:
        return f"distribution of {ref} clusters\ndataset {dataset}"


#calculate cluster distribution
def get_splits_merges(ref,method,dataset,organism,filter=None):
    #load data
    ref_df = pd.read_json(f"bench-{organism}/{dataset}-{ref}.json")
    method_df = pd.read_json(f"bench-{organism}/{dataset}-{method}.json")
    

    #filter data if provided
    if filter:
        # filter by lower and upper bounds (inclusive)
        # if only one value is specified use clusters up to hi or from lo on
        if isinstance(filter,tuple):
            lo = filter[0]
            hi = filter[1]
            if hi and lo:
                ref_df = ref_df[(ref_df["count"]>=lo)&(ref_df["count"]<=hi)]
            elif hi:
                ref_df = ref_df[ref_df["count"]<=hi]
            elif lo:
                ref_df = ref_df[ref_df["count"]>=lo]
        # filter by cluster size
        elif isinstance(filter,int):
            ref_df = ref_df[ref_df["count"]==filter]
        # filter by sccg or singleton (column specifier)
        elif isinstance(filter,str):
            ref_df = ref_df[ref_df[filter]=="yes"]

    #calculate hash for each cluster
    ref_df = ref_df.assign(hash=ref_df["locus"].apply(get_hash))
    method_df = method_df.assign(hash=method_df["locus"].apply(get_hash))

    #get identical clusters
    iden_hash = set(ref_df["hash"]).intersection(set(method_df["hash"]))
    identical = set(ref_df.query("hash in @iden_hash")["geneId"].tolist())

    #filter identical clusters
    ref_labels = ref_df.query("hash not in @iden_hash")[["locus","geneId","count"]]
    method_labels = method_df.query("hash not in @iden_hash")[["locus","geneId","count"]]
    #merge labels
    labels = ref_labels.explode("locus").merge(method_labels.explode("locus"),how="left",on="locus",suffixes=[ref,method])

    #get set of reference clusters that are split up
    split_id = labels.drop(["locus",f"count{ref}",f"count{method}"],axis=1).groupby(by=f"geneId{ref}").agg(lambda x: len(set_wo_nan(x))).reset_index()
    split_id = set(split_id[split_id[f"geneId{method}"]>1][f"geneId{ref}"].to_list())

    #get set of reference clusters that are merged
    merge_id = labels.drop("locus",axis=1).groupby(by=f"geneId{method}").agg(set = (f"geneId{ref}",set), len = (f"geneId{ref}",lambda x: len(set(x))), c_ref = (f"count{ref}","mean"), c_m = (f"count{method}","mean")).reset_index()
    
    if merge_id[(merge_id["len"]==1) & (merge_id["c_ref"]<merge_id["c_m"])]["set"].tolist():
        inflated = set.union(*merge_id[(merge_id["len"]==1) & (merge_id["c_ref"]<merge_id["c_m"])]["set"].tolist()).difference(split_id)
    else:
        inflated = set()

    if merge_id[(merge_id["len"]==1) & (merge_id["c_ref"]>merge_id["c_m"])]["set"].tolist():
        incomplete = set.union(*merge_id[(merge_id["len"]==1) & (merge_id["c_ref"]>merge_id["c_m"])]["set"].tolist()).difference(split_id)
    else:
        incomplete = set()

    if merge_id[merge_id["len"]>1]["set"].tolist():
        merge_id = set.union(*merge_id[merge_id["len"]>1]["set"].tolist())
    else:
        merge_id = set() 
    
    #get pure merges, pure splits and splits&merges
    splits = split_id.difference(merge_id)
    merges = merge_id.difference(split_id)
    both = split_id.intersection(merge_id)

    #get set of reference clusters that do not fall in any category (e.g. all locus tag not present in other method)
    other = set(ref_labels["geneId"]).difference((split_id.union(merge_id,inflated,incomplete)))

    return identical, splits, merges, both, inflated, incomplete, other


if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d","--dataset")
    parser.add_argument("-o","--organism")
    parser.add_argument("-r","--reference", help="reference method")
    parser.add_argument("-lo","--lower", help="filter by cluster size lower bound")
    parser.add_argument("-hi","--upper", help="filter by cluster size upper bound")
    parser.add_argument("-s","--size", help="filter by exact cluster size")
    parser.add_argument("-t","--type", help="filter by cluster type e.g. singletons or sccg")
    parser.add_argument("-p","--plot", action="store_true")
    parser.add_argument("-e","--export", action="store_true")

    args = parser.parse_args()

    dataset = args.dataset
    organism=args.organism
    ref = args.reference
    f = None
    fn = f"out/{ref}-{dataset}-all"
   

    print(f"reference:\t{ref}\ndataset:\t{dataset}")

    if args.lower and args.upper:
        f = (int(args.lower),int(args.upper))
        fn = f"out/{ref}-{dataset}-({args.lower},{args.upper})"
    elif args.lower:
        f = (int(args.lower),None)
        fn = f"out/{ref}-{dataset}-({args.lower},-)"
    elif args.upper:
        f = (None,int(args.upper))
        fn = f"out/{ref}-{dataset}-(-,{args.upper})"
    
    if args.size:
        f = int(args.size)
        fn = f"out/{ref}-{dataset}-{args.size}"
    if args.type:
        f = args.type
        fn = f"out/{ref}-{dataset}-{args.type}"

    if not os.path.exists("out"):
        os.makedirs("out")

    METHODS = ["roary","panacota","peppan","panaroo","panta","panX","panXX"]
    out_df = pd.DataFrame(columns=["method","identical","split","merge","split&merge","inflated","incomplete","other"])
    to_plot = pd.DataFrame(columns=["method","identical","split","merge","split&merge","inflated","incomplete","other"])

    for method in METHODS:
        if method==ref:
            continue
        print(f"comparing to {method}")
        try:
            i, s, m, b, infl, inco, o = get_splits_merges(ref,method,dataset,organism,f)
        except:
            print(f"Error ocurred: input file not available\n skipping {method}")
            continue
        out_df.loc[len(out_df)]=[method,list(i), list(s), list(m), list(b), list(infl), list(inco), list(o)]
        to_plot.loc[len(to_plot)]=[method,len(i),len(s),len(m),len(b),len(infl),len(inco),len(o)]
    
    #check if total cluster number is identical across methods
    n_clusters = set(to_plot.loc[:,to_plot.columns!="method"].sum(axis=1).tolist())
    if len(n_clusters)==1:
        total = n_clusters.pop()
    else:
        total = n_clusters
        print("inconsistent total size!")
    
    print("\nstatistics:")
    print(f"number of clusters: {total}")
    print(to_plot)
    
    if args.export:
        print(f"\nresults saved to: {fn}.json")
        out_df.to_json(f"{fn}.json",orient="records",index=False,indent=1)

    if args.plot:
        fig, ax = plt.subplots(1,1)
        to_plot.plot.bar(x="method",y=["split","merge","split&merge","inflated","incomplete","other"],stacked=True,ax=ax)

        ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left",title="types")
        ax.set_title(make_title(ref,dataset,f)+f" | {total} clusters")

        fig.tight_layout()
        print(f"\nplot saved to: {fn}.png")
        fig.savefig(f"{fn}.png", dpi=300)
