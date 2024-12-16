from BCBio import GFF
import json
import os
import argparse

def generate_prt(in_dir,out_dir):
    locus_dict = {}
    prt_dict = {}
    file_dict = {}
    n_file=1
    
    for file in os.listdir(in_dir):
        print(f"dealing with file {n_file}: {file}")
        n_file+=1
        strain = file[:-4]
        # create outdir
        if not os.path.exists(out_dir+"/prt"):
            os.makedirs(out_dir+"/prt")
        
        with open(in_dir+"/"+file,"r") as gff, open(out_dir+"/prt/"+file[:-3]+"prt","w") as prt:
            seq_id=None
            n_feat=0
            #filter by CDS and biological regions
            limit_info = dict(gff_type=["CDS","biological_region"])
            for rec in GFF.parse(gff,limit_info=limit_info):
                #maps record id to file name
                file_dict[rec.id]=strain
                #ensures prt entries all have the same species name
                if not seq_id:
                    seq_id = rec.id
                for feature in rec.features:
                    n_feat+=1
                    gff_name = feature.id
                    tag = feature.qualifiers.get("locus_tag")
                    if not tag:
                        # try to extract locus tag from linked parent
                        tag = feature.qualifiers.get("Parent")
                        if tag and ("." in tag[0]):
                            tag=[tag[0].split(".")[0]]
                    translation = feature.qualifiers.get("translation")

                    #if there is a locus tag add to locus dict
                    if gff_name and tag:
                        #maps locus tag to feature id
                        locus_dict[tag[0]]=gff_name

                    #if there is a translation available, write in prt file
                    if gff_name and translation:
                        prt_name = "_".join([seq_id,str(n_feat)])
                        prt.write(">%s\n%s\n" % (prt_name,translation[0]))
                        prt_dict[prt_name] = gff_name

    #export to json
    with open(out_dir+"/locus2id.json",'w') as out_file_1:
        json.dump(locus_dict, out_file_1, indent=1)

    with open(out_dir+"/prt2id.json",'w') as out_file_2:
        json.dump(prt_dict, out_file_2, indent=1)

    with open(out_dir+"/id2file.json",'w') as out_file_3:
        json.dump(file_dict, out_file_3, indent=1)

if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("in_dir")
    parser.add_argument("out_dir")
    args=parser.parse_args()

    if not args.out_dir:
        args.out_dir = args.in_dir + "/prt"
    
    generate_prt(args.in_dir,args.out_dir)
