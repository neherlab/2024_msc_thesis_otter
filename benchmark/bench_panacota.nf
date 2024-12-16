process to_lst{
    label 'q_mini'

    input:
    tuple val(dataset_name), path(prt)

    output:
    val dataset_name
    file("out/input.lst")
    path prt
    
    script:
    """
    mkdir -p out
    touch out/input.lst

    for file in ${prt}/*.prt; do
        echo \$(basename -s .prt \$file) >> out/input.lst
    done
    """
}

process run_panacota{
    label 'q_medium'

    conda "env/panacota.yml"

    input:
    val dataset_name
    file lst
    path prt

    output:
    tuple val(dataset_name), path("out/pangenome.lst.summary.txt"), path("out/pangenome.lst")

    script:
    """
    mkdir -p input
    cp ${prt}/*.prt input
    PanACoTA pangenome -l ${lst} -n sim -d input -o out -f pangenome.lst
    """

}

process format_panacota_output{
    label 'q_mini'

    conda "env/format.yml"

    input:
    tuple val(dataset_name), path(summary), path(pangenome), val(prt_lookup)    

    output:
    path("out/${dataset_name}-panacota.json")

    script:
    """
    mkdir -p out
    python3 $projectDir/scripts/format_panacota_output.py ${summary} ${pangenome} ${prt_lookup} ${dataset_name} out
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

    input_prt = Channel.fromPath(params.indir + "**/prt",type: 'dir').map{it -> [it.getParent().name, it]}

    id_lookup = Channel.fromPath(params.indir + "**/locus2id.json",type: 'file').map{it -> [it.getParent().name, it]}

    prt_lookup = Channel.fromPath(params.indir + "**/prt2id.json",type: 'file').map{it -> [it.getParent().name, it]}
    
    input_prt | to_lst | run_panacota
    run_panacota.out.join(prt_lookup) | format_panacota_output  

    out_ch = Channel.fromPath(params.outdir,type: 'dir')

    format_panacota_output.out.combine(out_ch) | export
}