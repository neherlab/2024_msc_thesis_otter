process run_panX{
    label 'q_medium_16h'

    conda "env/panX.yml"

    input: 
    tuple val(species_name), path(gbk)

    output:
    tuple val(species_name), path("out/vis/geneCluster.json")


    script:
    """
    mkdir -p out
    cp ${gbk}/*.gbk out
    python $projectDir/panX/panX.py -fn out -sl ${species_name} -t 8 -dmdc
    """
}

process format_panX_output{
    label 'q_mini'

    conda "env/format.yml"

    input:
    tuple val(dataset_name), path(input_json), path(id_lookup)

    output:
    path("out/${dataset_name}-panX.json")

    script:
    """
    mkdir -p out
    python3 $projectDir/scripts/format_panX_output.py ${input_json} ${dataset_name} ${id_lookup} out
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
    // multiple datasets
    // params.indir = "/scicore/home/neher/GROUP/data/2024_panXX_bench/Ecoli-200"
    // params.outdir = "bench-Ecoli-n/"

    // input_gbk = Channel.fromPath(params.indir + "**/gbk",type: 'dir').map{it -> [it.getParent().name, it]}

    // id_lookup = Channel.fromPath(params.indir + "**/locus2id.json",type: 'file').map{it -> [it.getParent().name, it]}

    // prt_lookup = Channel.fromPath(params.indir + "**/prt2id.json",type: 'file').map{it -> [it.getParent().name, it]}

    
    // single dataset
    params.indir = "/scicore/home/neher/GROUP/data/2024_panXX_bench/Ecoli-500/500_1"
    params.outdir = "bench-Ecoli-s/"

    input_gbk = Channel.fromPath(params.indir + "/gbk",type: 'dir').map{it -> [it.getParent().name, it]}

    id_lookup = Channel.fromPath(params.indir + "/locus2id.json",type: 'file').map{it -> [it.getParent().name, it]}

    input_gbk | run_panX 
    run_panX.out.join(id_lookup) | format_panX_output

    out_ch = Channel.fromPath(params.outdir,type: 'dir')

    format_panX_output.out.combine(out_ch) | export
}