# Comparative evaluation of bacterial pan-genome analysis tools
## Master's Thesis of Benjamin Otter

# Benchmark

### input data
All input data is located at `/scicore/home/neher/GROUP/data/2024_panXX_bench/`
The input directory is set in the individual workflow files:
```
params.indir = "/scicore/home/neher/GROUP/data/2024_panXX_bench/Ecoli-500"
```

### output
The workflows automatically produce formatted output files. 
The output directory can be set in the individual workflow files:
```
params.outdir = "bench-Ecoli/"
```

### config
`nextflow.config` contains nextflow specific configuration. Adjustments to conda environments, docker images, executors as well as CPU and memory settings for each process can be made there.

### running the methods on the cluster
```
nextflow run bench_method.nf -profile cluster -with-report -with-trace
```

### panXX
PanXX is run seperately:
```
nextflow run panXX.nf -profile cluster -with-report spneumo.html -with-trace --batchsize=500 --in=/scicore/home/neher/GROUP/data/2024_panXX_bench/Spneumo/gbk --out=batch/spneumo --run="Spneumo"
```

Output files can be formatted using `format_panXX.nf`. The required parameters are:
- `params.indir`: path to panXX output files
- `params.ref`: path to original input dataset
- `params.outdir`: output directory for formatted file

The parameters can be set directly in the workflow file
```
params.indir = "/panXX/batch/spneumo"
params.ref = "/scicore/home/neher/GROUP/data/2024_panXX_bench/Spneumo"
params.outdir = "bench-Spneumo/"
```


