#!/bin/bash

for n in 100 200 300 400 500; do
  for i in 1 2 3 4 5; do

    echo "working on ${n}_${i}"

    X=/scicore/home/neher/GROUP/data/2024_panXX_bench/Escherichia_coli/gbk
    # Y=/scicore/home/neher/GROUP/data/2024_panXX_bench/Escherichia_coli/prt
    Z=/scicore/home/neher/GROUP/data/2024_panXX_bench/Escherichia_coli/gff

    A=/scicore/home/neher/GROUP/data/2024_panXX_bench/Ecoli-${n}/${n}_${i}/gbk
    # B=/scicore/home/neher/GROUP/data/2024_panXX_bench/Ecoli-1000/1000_1/prt
    C=/scicore/home/neher/GROUP/data/2024_panXX_bench/Ecoli-${n}/${n}_${i}/gff

    mkdir -p $A $C

    mapfile -t files < <(find "$X" -type f -name "*.gbk" -print | shuf -n $n -)

    for f in "${files[@]}"; do
      cp "$f" "$A"
      bn=$(basename "$f" ".gbk")
      # cp "$Y/$bn.prt" "$B"
      cp "$Z/$bn.gff" "$C"
    done
  done
done
