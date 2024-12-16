process run_roary{
    label 'q_medium_6h'

    conda "env/roary.yml"

    input: 
    tuple val(dataset_name), path(gff)

    output:
    tuple val(dataset_name), path("out/gene_presence_absence.csv")

    script:
    """
    roary -p 8 -f out ${gff}/*.gff
    """

}

process format_roary_output{
    label 'q_mini'

    conda "env/format.yml"

    input:
    tuple val(dataset_name), path(input_csv)

    output:
    path("out/${dataset_name}-roary.json")

    script:
    """
    mkdir -p out
    python3 $projectDir/scripts/format_roary_output.py ${input_csv} ${dataset_name} out
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

    input_gff | run_roary | format_roary_output

    out_ch = Channel.fromPath(params.outdir,type: 'dir')

    format_roary_output.out.combine(out_ch) | export
}