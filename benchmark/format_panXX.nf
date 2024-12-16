process format_panXX{
    label 'q_mini'

    conda "env/format.yml"

    input:
    tuple val(dataset_name), path(input_json), path(id_lookup)

    output:
    path("out/${dataset_name}-panXX.json")

    script:
    """
    mkdir -p out
    python3 $projectDir/scripts/format_panXX_output.py ${input_json} ${dataset_name} ${id_lookup} out
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
    /**
     * for multiple datasets 
     */

    // params.indir = "/scicore/home/neher/ottben00/Projects/panXX/batch/"
    // params.ref = "/scicore/home/neher/GROUP/data/2024_panXX_bench/"
    // params.outdir = "bench-SimPan/"

    // input_json = Channel.fromPath(params.indir + "*/geneCluster.json",type: 'file').map{it -> [it.getParent().name, it]}

    // id_lookup = Channel.fromPath(params.ref + "*/*/locus2id.json",type: 'file').map{it -> [it.getParent().name, it]}
    
    // input_json.join(id_lookup) | format_panXX

    // out_ch = Channel.fromPath(params.outdir,type: 'dir')

    // format_panXX.out.combine(out_ch) | export



    /**
     * for single dataset [prochlorococcus, spneumo or campy]
     */
    
    params.indir = "/scicore/home/neher/ottben00/Projects/panXX/batch/prochlorococcus"
    params.ref = "/scicore/home/neher/GROUP/data/2024_panXX_bench/Prochlorococcus"
    params.outdir = "bench-Prochlorococcus/"

    input_json = Channel.fromPath(params.indir + "/*/geneCluster.json",type: 'file').map{it -> [it.getParent().name, it]}

    id_lookup = Channel.fromPath(params.ref + "/locus2id.json",type: 'file').map{it -> [it.getParent().name, it]}

    input_json.join(id_lookup) | format_panXX

    out_ch = Channel.fromPath(params.outdir,type: 'dir')

    format_panXX.out.combine(out_ch) | export

}