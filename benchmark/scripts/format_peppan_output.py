import pandas as pd
import click
import json

@click.command()
@click.argument("input_csv", type=click.File("r"))
@click.argument("input_gff", type=click.File("r"))
@click.argument("id_lookup", type=click.File("r"))
@click.argument("dataset_name")
@click.argument("dst")
def format(input_csv,input_gff,id_lookup,dataset_name,dst):
    all_data = pd.read_csv(input_csv)
    max_count = all_data.shape[1]-1
    string_dict = pd.read_csv(input_gff,delimiter="\t",skiprows=2,header=None).iloc[:,8]
    locus_dict={}

    lookup = json.load(id_lookup)

    for line in string_dict:
        fields = line.split(";")
        peppan_name = None
        tag = None

        for field in fields:
            if field.startswith("ID="):
                #PEPPAN_g_XXXXX
                peppan_name = field[3:]
            elif field.startswith("old_locus_tag="):
                tag = field.split(":")[0][14:]
        #if we have both variables, then the locus tag is in lookup table
        if peppan_name and tag:
            #if a gene has no locus tag, peppan assigns the id to this field, which is needed in the final output, continue to the next line
            if "." in tag:
                print(f"{tag} is already formatted.")
                locus_dict[peppan_name]=tag
                continue
            # for the case if some locus tags unexpectedly went missing, print warning but add 'none' anyway
            elif not lookup.get(tag):
                print(f"{tag} should be in idlookup json")
            locus_dict[peppan_name]=lookup.get(tag)
        else:
            #peppan added this gene, no way to assign a id to compare with other methods
            locus_dict[peppan_name]=peppan_name

    locuslist=all_data.iloc[:,1:].apply(combine_columns,locus_dict=locus_dict,max_count=max_count,axis=1)

    metadata=pd.DataFrame(locuslist.tolist(),columns=["locus","count","dupli","sccg","singleton"])
    metadata["geneId"]=metadata.index+1

    metadata=metadata[["geneId","count","dupli","sccg","singleton","locus"]]

    metadata.to_json(f"{dst}/{dataset_name}-peppan.json",orient='records',indent=1)



def combine_columns(row,locus_dict,max_count):
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
                    locus_list.append(locus_dict[t])
                    l+=1
                dupli="yes"
            else:
                l+=1
                locus_list.append(locus_dict[tag])
    if s==1:
        singleton="yes"
    if (l==max_count) and (dupli=="no"):
        sccg="yes"

    return locus_list, l, dupli, sccg, singleton

format()