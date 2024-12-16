process export{
    label 'q_mini'

    input:
    tuple val(key), val(format), path(filelist), val(out)

    output:

    script:
    """
    mkdir -p ${out}
    outdir=${out}/${key}
    mkdir -p \$outdir/gbk \$outdir/gff
    
    for file in ${filelist}
    do  
        ext=\${file##*.}
        if [ "\$ext" != "fna" ] && [ "\$ext" != "${format}" ]; then
            cp -n \$file \$outdir/\$ext
        fi
    done
    """
}

def determineFormatClass(List<String> fileList) {
    def uniqueFormats = fileList.toSet().toSorted()
    return uniqueFormats.join('+')
}

process gbk2prt{
    label 'q_mini_long'

    conda 'env/biopython-click.yml'

    input:
    tuple val(key), val(format), path(files)

    output:
    tuple val(key), val(format), path("out/*.prt"), emit: new_file

    """
    mkdir -p out

    for file in ${files}
    do
        filename=\$(basename -s .gbk \$file)
        touch out/\$filename.prt
        python3 $projectDir/scripts/gbk2prt.py \$filename.gbk out/\$filename.prt
    done   
    """
}

process gff2gbk{
    label 'q_mini'

    conda 'env/bioperl-emboss.yml'

    input:
    tuple val(key), val(format), path(files)

    output:
    tuple val(key), val(format), path("out/*.gbk"), emit: new_file

    script:
    """
    mkdir -p out intermediate

    for file in ${files}
    do
        if [[  \$file == *.gff ]]; then

            filename=\$(basename -s .gff \$file) 
            fasta=\$filename.fna

            seqret -sequence \$fasta -feature -fformat gff -fopenfile \$file -osformat genbank -osname_outseq \$filename -osdirectory_outseq intermediate -auto
            $projectDir/scripts/gb-add-trans intermediate/\$filename.genbank > out/\$filename.gbk
        fi
    done  
    """
}

process gbk2gff_fna{
    label 'q_mini'

    conda 'env/bioperl-emboss.yml'

    input:
    tuple val(key), val(format), path(files)

    output:
    tuple val(key), val(format), path("out/*.{gff,fna}"), emit: new_file

    script:
    """
    mkdir -p out
    touch header.txt
    echo "##FASTA" > header.txt
    for file in ${files}
    do
        filename=\$(basename -s .gbk \$file) 
        seqret -sformat genbank -sequence \$file -feature -offormat gff3 -osformat fasta -osname_outseq \$filename -osdirectory_outseq out -ofname_outseq \$filename.gff -ofdirectory_outseq out -supper_sequence -auto 
    done 

    for file in out/*.fasta 
    do
        filename=\$(basename -s .fasta \$file)
        cat header.txt out/\$filename.fasta >> out/\$filename.gff
        mv -- \$file out/\$filename.fna
    done

    """
}

process gbk2gff3{
    label 'q_mini'

    conda 'env/biocode.yml'

    input:
    tuple val(key), val(format), path(files)

    output:
    tuple val(key), val(format), path("out/*.gff"), emit: new_file

    script:
    """
    mkdir -p out
    
    for file in ${files}
    do
        filename=\$(basename -s .gbk \$file) 
        $projectDir/scripts/convert_genbank_to_gff3.py -i \$file -o out/\$filename.gff --with_fasta
    done 
    """
}

process bp_gbk2gff3{
    label 'q_mini_long'

    input:
    tuple val(key), val(format), path(files)

    output:
    tuple val(key), val(format), path("out/*.gff"), emit: new_file

    script:
    """
    mkdir -p out
    $projectDir/scripts/bp_genbank2gff3 --out out ${files}

    for file in out/*.gff
    do
        filename=\$(basename -s .gff \$file)
        strain=\$(basename -s .gbk \$filename)
        mv -- \$file out/\$strain.gff
    done
    """
}


workflow {

    params.in = "/scicore/home/neher/GROUP/data/2024_panXX_bench/SimPan-15"
    params.out = "/scicore/home/neher/GROUP/data/2024_panXX_bench/SimPan-15"

    input_ch = Channel.fromPath(params.in + '**/*.{gff,gbk,fna}')

    mapped_ch = input_ch.map { file ->
        def pathString = file.toString()
        def parts = pathString.split('/')
        def species = parts[-3]
        def format = file.name.tokenize('.').last()
        return [format, species, file]
    }

    mapped_ch.groupTuple(by: [1])
        .map{formatlist, species, files -> [species, determineFormatClass(formatlist), files]}
        .branch {
            gbk: it[1] == "gbk"
            gff: it[1] == "gff"
            fna_gff: it[1] == "fna+gff"
            }
        .set{sorted_ch}

    sorted_ch.fna_gff | gff2gbk

    // conversion to gff3
    sorted_ch.gbk | gbk2gff_fna
    gbk2gff_fna.out.mix(gff2gbk.out, sorted_ch.fna_gff, sorted_ch.gbk).groupTuple(by:[0,1]).set{export_ch}

    // alternative conversion to gff3 (used for S.pneumniae dataset)
    // sorted_ch.gbk | bp_gbk2gff3
    // bp_gbk2gff3.out.mix(gff2gbk.out, sorted_ch.fna_gff, sorted_ch.gbk).groupTuple(by:[0,1]).set{export_ch}
    

    // export gff and gbk files
    export_ch.map{species,format,filelist -> [species, format, filelist.flatten(), params.out]} | export

}