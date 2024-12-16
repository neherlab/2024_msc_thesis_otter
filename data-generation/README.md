# Comparative evaluation of bacterial pan-genome analysis tools
## Master's Thesis of Benjamin Otter

# Data Generation

### general workflow
- start with a set of .gbk or .gff files, or generate a simulated dataset
- run `main.nf` to generate gff, and gbk files (works on multiple species folders)
- optional subsampling individual runs per species with `scripts/subsample.sh`
- for each (subsampled) species run `write_prt.nf` to generate prt files and lookup tables

### config
`nextflow.config` contains nextflow specific configuration. Adjustments to conda environments, docker images, executors as well as CPU and memory settings for each process can be made there.

### generating a simulated dataset
Simulated datasets are created with `sim_data.nf`.
Set up parameters:
```
params.out = "/scicore/home/neher/GROUP/data/2024_panXX_bench/SimPan-50"

idenOrtho = Channel.value(0.7)
idenPara = Channel.value(0.5)
prefix = Channel.value("70_50")
n_genomes = Channel.value(50)
```
Run the workflow:
```
nextflow run sim_data.nf -profile cluster
```

### prepare a dataset: format conversion
GenBank and GFF files are generated using `main.nf`.
Set up parameters:
```
params.in = "/scicore/home/neher/GROUP/data/2024_panXX_bench/SimPan-50"
params.out = "/scicore/home/neher/GROUP/data/2024_panXX_bench/SimPan-50"
```
Run the workflow:
```
nextflow run main.nf -profile cluster
```

### prepare a dataset: prt files and lookup tables
prt files and lookup tables are generated using `write_prt.nf`.
Set up parameters:
```
params.in = "/scicore/home/neher/GROUP/data/2024_panXX_bench/SimPan-50/70_50"
```

Run the workflow:
```
nextflow run write_prt.nf -profile cluster
```

The following lookup tables are created:
- `locus2id.json` maps locus tags to id, essential for panX, panXX and PEPPAN formatting
- `prt2id.json` maps protein names to id, essential for PanACoTA formatting
- `id2file.json` maps contig names to original file name, essential for retrieving sequences in `cluster_to_tree.py`



