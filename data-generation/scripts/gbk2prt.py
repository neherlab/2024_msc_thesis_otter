from Bio import SeqIO
import click
import json
import os

@click.command()
@click.argument("input_gbk", type=click.File("r"))
@click.argument("output_prt", type=click.File("w"))
@click.argument("prt2id_json")
def gbk2prt(input_gbk, output_prt, prt2id_json):
    for seq_record in SeqIO.parse(input_gbk, "genbank"):      
        for seq_feature in seq_record.features:
            if seq_feature.type=="CDS":
                #extract relevant feautures
                translation = seq_feature.qualifiers.get('translation')
                tag = seq_feature.qualifiers.get('locus_tag')
                #skip genes that have no translation or locus_tag
                if translation and tag:
                    #write prt file
                    output_prt.write(">%s\n%s\n" % (tag,translation[0]))
    
gbk2prt()