## SCiMS-paper

This repository contains reproducible scripts for the SCiMS manuscript.

The associated preprint will be added to the repository once it is published.

<!-- TOC start -->

### Contents:
- [Environment setup](#environment-setup)
- [Part 1: Simulation](#part-1-simulation)
- [Part 2: Human Microbiome Project (HMP) dataset](#part-2-human-microbiome-project-hmp-dataset)
- [Part 3: Mouse Metagenomic Dataset](#part-3-mouse-metagenomic-dataset)
- [Part 4: Chicken Metagenomic Dataset](#part-4-chicken-metagenomic-dataset)
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
snakemake -s 01_simulation/01_simulation_create_reads.smk --cores 4 --use-conda
```
Expected outputs:
- Simulated reference genomes (male and female) in `./data/simulated_ref`
- Simulated FASTQ reads in `./data/simulated_reads`

2. Map simulated reads and process BAM files:

```bash
snakemake -s 01_simulation/02_simulation_map_and_process_simulated_reads.smk --cores 4 --use-conda
```
Expected outputs:
- Sorted, indexed, and duplicate-removed BAM files (sorted.rmdup.bam) in `./data/mapped_reads`

3. Downsample the BAM files:
```bash
snakemake -s 01_simulation/03_simulation_downsample_reads.smk --cores 4 --use-conda
```
Expected outputs:
- Downsampled BAM files (*.1000x.bam) in `./data/mapped_reads`
- Index stats files (*.1000x.idxstats) in `./data/mapped_reads`

4. Run SCiMS on downsampled simulated data:
```bash
bash 01_simulation/04_simulation_scims.sh
```
Expected outputs:
- SCiMS results in `./results/scims`

5. Run Rx, Ry, and BeXY on downsampled simulated data:
```bash
bash 01_simulation/05_simulation_rxry.sh
bash 01_simulation/06_simulation_bexy.sh
```

7. Generate the Figure 2:
```bash
python3 ./scripts/figure2.py
```

Expected outputs:
- Figures in `./figures`




