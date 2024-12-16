process gff2prt{
    label 'q_mini_long'

    conda 'env/bcbio.yml'

    input:
    path(in_dir)

    output:
    path("out/prt/")
    path("out/")

    """
    mkdir -p out
    python3 $projectDir/scripts/gff2prt.py ${in_dir} out
    """
}

process export_prt{
    label 'q_mini'

    input:
    tuple val(out), val(prt), val(json)

    output:

    script:
    """
    mkdir -p ${out}/prt
    
    for file in ${prt}/*.prt
    do  
        cp \$file ${out}/prt
    done

    for file in ${json}/*.json
    do  
        cp \$file ${out}
    done
    """
}

workflow {

    params.in = "/scicore/home/neher/GROUP/data/2024_panXX_bench/SimPan-15/90_80"

    input_ch = Channel.fromPath(params.in + "/gff", type: 'dir')
    out_ch =  Channel.fromPath(params.in, type: 'dir')

    input_ch | gff2prt
    
    out_ch.concat(gff2prt.out).toList() | export_prt
}

