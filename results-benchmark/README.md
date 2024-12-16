# Comparative evaluation of bacterial pan-genome analysis tools
## Master's Thesis of Benjamin Otter

# results benchmark
### raw data
`trace` contains the trace files per run and tool. This data was used to compute peak memory usage.
`ecoli-cputime.txt` contains the reported CPU-times per run and tool
`panXX-trace-500-re-nf.txt` contains the results for one single panXX re-run on a 500 ecoli dataset. This data was used to assess the time and memory usage of individual panXX processes

### generate plots
Plots are generated with `plot_cpu_hours.ipynb`

### additional scripts
`pyvenn/create_venn5.ipynb` is used to create a venn-5 diagram to depict the pan-genome-concept
`panXX-trace.ipynb` is used to assess the time and memory usage of each individual panXX process
`pan-genome-simulation` is used to plot a theoretical core- and pan-genome-size for increasing number of genomes 
