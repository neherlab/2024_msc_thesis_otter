import argparse
import pandas as pd
import distinctipy
import json
import subprocess
import os

from BCBio import GFF

from phytreeviz import TreeViz
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle


# --------------------------- #
#          methods
# --------------------------- #

#extract protein sequence from gff file
def extract_sequence(tag,dataset):
    datapath = f"/scicore/home/neher/GROUP/data/2024_panXX_bench/{dataset}"
    with open(f"{datapath}/id2file.json") as fh:
        id2file = json.load(fh)
    strain = id2file.get(tag.split(".")[0])
    infile = f"{datapath}/gff/{strain}.gff"
    limit_info = dict(gff_type=["databank_entry","CDS","biological_region"])
    with open(infile) as fh:
        for rec in GFF.parse(fh,limit_info=limit_info):
            if rec.features:
                for feat in rec.features:
                    if feat.type == "databank_entry":
                        meta = feat.qualifiers.get("organism")[0]
                    elif feat.id == tag:
                        if feat.qualifiers.get("translation"):
                            return [feat.qualifiers.get("translation")[0],meta]
                        else:
                            return ["not available",meta]

#extract geneId given a locus tag                
def get_label(tag,df):
    try:
        return df[df["locus"]==tag].geneId.values[0]
    except: 
        return "not found"

def get_species(meta):
    if "coli" in meta:
        return "C.coli"
    elif "jejuni" in meta:
        return "C.jejuni"
    else:
        print("organism info not available")
        print(meta)
        return None

#mafft
def align_sequences(input_file, output_file):
    mafft_command = f"mafft --auto {input_file} > {output_file}"
    subprocess.call(mafft_command, shell=True)

#FastTree
def build_tree(input_file, output_file):
    fasttree_command = f"FastTree {input_file} > {output_file}"
    subprocess.call(fasttree_command, shell=True)


if __name__=="__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-d","--dataset")
    parser.add_argument("-o","--organism")
    parser.add_argument("-i","--geneId", type=int, help="geneId of reference cluster")
    parser.add_argument("-r","--reference", help="reference method")
    parser.add_argument("-l","--label", help="label method")

    args = parser.parse_args()

    dataset = args.dataset
    organism = args.organism
    geneId = args.geneId
    ref = args.reference
    method = args.label

    if not os.path.exists("out"):
        os.makedirs("out")
    
    if not os.path.exists("out/tmp"):
        os.makedirs("out/tmp")
    

    # ------------------------------------------ #
    #   get list of locus tags based on geneId
    #   and export results to fasta file
    # ------------------------------------------ #
    ref_df = pd.read_json(f"bench-{organism}/{dataset}-{ref}.json")
    method_df = pd.read_json(f"bench-{organism}/{dataset}-{method}.json")

    #select cluster of interest
    cluster = ref_df[ref_df["geneId"]==geneId]
    sccg = cluster["sccg"].values[0]
    sing = cluster["singleton"].values[0]
    dupli = cluster["dupli"].values[0]
    cluster = cluster[["geneId","locus"]].explode("locus")
    size = len(cluster)
    print(f"Cluster of interest\ngeneId:\t {geneId}\nmethod:\t {ref}\ndataset: {dataset}\ncount:\t {size}\n")

    #get sequences for each tag
    print("retrieving sequences")
    if organism=="Ecoli":
        dataset_path=f"E-coli/{dataset[:-2]}/{dataset}"
    else:
        dataset_path = dataset
    
    cluster[["sequence","metadata"]] = cluster["locus"].apply(extract_sequence,dataset=dataset_path).to_list()
    
    if organism == "Campy":
        cluster["metadata"] = cluster["metadata"].apply(get_species)

    #filter out any entries without a sequence
    to_align = cluster[~(cluster["sequence"]=="not available")]
    extra = cluster[cluster["sequence"]=="not available"]

    #get method label for annotation
    print(f"retrieving {method} labels")
    method_labels = method_df[["geneId","locus"]].explode("locus")
    to_align = to_align.assign(label=to_align["locus"].apply(get_label,df=method_labels))

    #export to fasta file
    outfile = "out/tmp/align_in.fasta"
    with open(outfile,"w") as f:
        for _,tag,seq,label in to_align[["locus","sequence","label"]].itertuples():
            f.write(f">{tag}\n{seq}\n")
    print(f"exported results to {outfile}\n")

    # ----------------------------- #
    #       align with mafft
    # ----------------------------- #
    infile = "out/tmp/align_in.fasta"
    outfile = "out/tmp/align_out.fasta"
    with open(infile) as fh:
        #skip header
        next(fh)
        #get len of first sequence (-1 for trailing newline)
        seq_len = len(next(fh))-1
    align_sequences(infile, outfile)
    print("Alignment done!\n")

    # ----------------------------- #
    #          build tree
    # ----------------------------- #
    infile = "out/tmp/align_out.fasta"
    outfile = "out/tmp/cluster_tree.nwk"
    build_tree(infile, outfile)
    print("\nBuilt tree!\n")

    # ----------------------------- #
    #        visualize tree
    # ----------------------------- #
    infile = "out/tmp/cluster_tree.nwk"
    outfile = f"out/{organism}_{ref}_{geneId}_{method}.png"
    height = 0.2
    tv = TreeViz(infile,height=height,align_leaf_label=True,leaf_label_size=10)
    tv.show_scale_bar()

    labels = to_align[["label","locus"]].groupby(by="label").agg(list).reset_index()
    colors = distinctipy.get_colors(len(labels),pastel_factor=0.2)
    handles = []
    for i,row in labels.iterrows():
        c = colors[i]
        label = row["label"]
        if label == "not found":
            c="black"
            handles.append(Patch(label=label, color=c))
        else:
            count = method_df[method_df["geneId"]==label]["count"].values[0]
            handles.append(Patch(label=f"{label} [{count}]", color=c))
            
        for l in row["locus"]:
            tv.marker(l, marker="D", color=c, descendent=False)

    #add missing genes
    if len(extra)>0:
        handles.append(Rectangle((0,0),1,1, fill=False, edgecolor=None, visible=False, label="\nno sequence"))
        handles.append(Rectangle((0,0),1,1, fill=False, edgecolor=None, visible=False, label="\n".join(extra["locus"])))


    i=0
    group_colors = ["darkslategrey","darkred","midnightblue"]
    handles.append(Rectangle((0,0),1,1, fill=False, edgecolor=None, visible=False, label="\norganism"))
    for _,group in to_align.groupby(by="metadata"):
        handles.append(Patch(label=_, color=group_colors[i%10]))
        for l in group.locus.to_list():
            tv.set_node_label_props(l,color=group_colors[i%10])
        i+=1

    fig,ax = plt.subplots(figsize=(8,max(height*len(to_align),5)))
    ax.set_title(f"cluster tree {organism} | {ref} geneId {geneId} | dataset {dataset}\n count {size} | singleton {sing} | sccg {sccg} | duplication {dupli} | protein length {seq_len}")

    tv.plotfig(ax=ax)
    fig.legend(
        handles=handles,
        title=f"{method} geneId",
        frameon=True,
        bbox_to_anchor=(1,0.985),
        loc="upper left"
    )

    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    print(f"Figure saved to: {outfile}\n")