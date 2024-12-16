process run_peppan{
    label 'q_medium_16h'

    conda "env/peppan.yml"

    input: 
    tuple val(dataset_name), path(gff)

    output:
    tuple val(dataset_name), path("PEPPAN.PEPPAN.gene_content.csv"), path("PEPPAN.PEPPAN.gff")

    script:
    """
    PEPPAN ${gff}/*.gff
    PEPPAN_parser -g PEPPAN.PEPPAN.gff
    """

}

process format_peppan_output{
    label 'q_mini'

    conda "env/format.yml"

    input:
    tuple val(dataset_name), path(input_csv), path(input_gff), path(id_lookup)

    output:
    path("out/${dataset_name}-peppan.json")

    script:
    """
    mkdir -p out
    python3 $projectDir/scripts/format_peppan_output.py ${input_csv} ${input_gff} ${id_lookup} ${dataset_name} out
    """
}

process export{
    label "q_mini"

    input:
    tuple path(file_list), path(dst) 

    output:
    

    script:
    """
    mkdir -p $projectDir/${dst}
    for file in ${file_list}
    do  
        cp \$file $projectDir/${dst}
    done
    """
}


workflow{
    params.indir = "/scicore/home/neher/GROUP/data/2024_panXX_bench/Ecoli-500"
    params.outdir = "bench-Ecoli-n/"

    input_gff = Channel.fromPath(params.indir + "**/gff",type: 'dir').map{it -> [it.getParent().name, it]}

    id_lookup = Channel.fromPath(params.indir + "**/locus2id.json",type: 'file').map{it -> [it.getParent().name, it]}

    prt_lookup = Channel.fromPath(params.indir + "**/prt2id.json",type: 'file').map{it -> [it.getParent().name, it]}
    
    input_gff | run_peppan 
    run_peppan.out.join(id_lookup) | format_peppan_output

    out_ch = Channel.fromPath(params.outdir,type: 'dir')

    format_peppan_output.out.combine(out_ch) | export
}