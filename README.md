## SCiMS-paper

This repository contains reproducible scripts for the SCiMS manuscript.

The associated preprint will be added to the repository once it is published.

<!-- TOC start -->

### Contents:
- [Environment setup](#environment-setup)
- [Simulation](#simulation)
- [Human Microbiome Project (HMP) dataset](#human-microbiome-project-hmp-dataset)
- [Mouse Metagenomic Dataset](#mouse-metagenomic-dataset)
- [Chicken Metagenomic Dataset](#chicken-metagenomic-dataset)
- [References](#references)

<!-- TOC end -->

### Environment setup
1. Clone the repository
```bash
git clone https://github.com/hanhntran/SCiMS-paper.git
cd SCiMS-paper
```

2. Install mamba: if you haven't installed mamba yet, use the following command to install it:
```bash
conda install -y -c conda-forge -c bioconda mamba
```

3. Create and activate the conda environment:
```bash
mamba env create -n scims-env -f ./envs/environment.yaml
mamba activate scims-env
```

### Part 1: Simulation
1. Generate the simulation data:
```bash
cd 01_simulation
snakemake -s 01_simulation_create_reads.smk --cores 4 --use-conda
```
Expected outputs:
- Simulated reference genomes (male and female) in `./data/simulated_ref`
- Simulated FASTQ reads in `./data/simulated_reads`

2. Map simulated reads and process BAM files:

```bash
snakemake -s 02_simulation_map_and_process_simulated_reads.smk --cores 4 --use-conda
```
Expected outputs:
- Sorted, indexed, and duplicate-removed BAM files (sorted.rmdup.bam) in `./data/mapped_reads`

3. Downsample the BAM files:
```bash
snakemake -s 03_simulation_downsample_reads.smk --cores 4 --use-conda
```
Expected outputs:
- Downsampled BAM files (*.1000x.bam) in `./data/mapped_reads`
- Index stats files (*.1000x.idxstats) in `./data/mapped_reads`

4. Run SCiMS on downsampled simulated data:
```bash
bash 04_simulation_scims.sh
```
Expected outputs:
- SCiMS results in `./results/scims`






