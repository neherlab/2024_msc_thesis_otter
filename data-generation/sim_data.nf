process run_simpan{
    label 'q_medium_long'

    input:
    tuple val(idenOrtho), val(idenPara), val(n_genomes), val(prefix)

    output:
    tuple val(prefix), val(n_genomes), path("*.fna"), path("*.gff"), path("${prefix}_*.tbl")

    script:
    """
    sb --genomeNum ${n_genomes} --aveSize 4500 --nBackbone 4000 --nCore 3500 --nMobile 10000 --rec 0.1 --idenOrtholog ${idenOrtho} --idenParalog ${idenPara} --prefix ${prefix}
    """
}

process make_reference{
    label 'q_mini'

    conda 'env/pandas.yml'

    input:
    tuple val(prefix), val(n_genomes), path(fna), path(gff), path(tbl)

    output:
    val(prefix)
    path("out/${prefix}-reference.json")
    path(fna)
    path(gff)
    path(tbl)

    script:
    """
    mkdir -p tmp
    # strip and combine files
    for filename in *.tbl; do
        awk '(NR+3) % 6 == 0 || (NR+2) % 6 == 0' \$filename | awk '{ print \$NF }' | sed '\$!N;s/\\n/\\t/' > tmp/s-\$filename
        cat tmp/s-\$filename >> combined.tbl
    done

    mkdir -p out
    # python make_reference prefix combined.tbl path to tmp files gene.tsv
    python $projectDir/scripts/create_reference.py ${n_genomes} ${prefix} combined.tbl tmp out
    """

}


process export{
    label 'q_mini'

    input:
    tuple path(out), val(prefix), path(ref_json), path(fna), path(gff), path(tbl)

    output:

    script:
    """
    mkdir -p ${out}/${prefix}/{fna,gff,tbl}
     
    for file in ${fna}
    do  
        cp \$file ${out}/${prefix}/fna
    done 

    for file in ${gff}
    do  
        cp \$file ${out}/${prefix}/gff
    done

    for file in ${tbl}
    do  
        cp \$file ${out}/${prefix}/tbl
    done

    cp ${ref_json} ${out}/${prefix}
    """
}


workflow{
    params.out = "/scicore/home/neher/GROUP/data/2024_panXX_bench/SimPan-15"
    out_dir = Channel.fromPath(params.out, type: 'dir')
    
    idenOrtho = Channel.value(0.7)
    idenPara = Channel.value(0.5)

    prefix = Channel.value("70_50")
    n_genomes = Channel.value(50)

    idenOrtho.concat(idenPara, n_genomes, prefix).set{input_ch}
    
    input_ch.toList() | run_simpan | make_reference

    out_dir.concat(make_reference.out).toList() | export

}