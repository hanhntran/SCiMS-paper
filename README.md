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
> [!WARNING] 
> The downsampling step will generate ~72,000 files and will take several days to complete and memory-intensive (~300GB). \
> Alternatively, you can run this pipeline in parallel using a Snakemake plugin 'snakemake-executor-plugin-slurm`. \
> For example, to run the pipeline in parallel, you can use the following commands. With these settings, the pipeline will run in parallel on ~25 jobs, and each job will take 3-4 hours to complete.

```bash
pip install snakemake-executor-plugin-slurm

snakemake -s 01_simulation/03_simulation_downsample_reads.smk \
          --executor cluster-generic \
          --cluster-generic-submit-cmd "sbatch --partition=open --account=open --time=12:00:00 --nodes=1 --ntasks=1 --mem=20GB" \
          --jobs 60 \
          --groups downsample_bam=group_03 index_downsampled_bam=group_03 generate_idxstats=group_03 \
          --group-components group_03=60 \
          --rerun-incomplete \
          --latency-wait 60 
```

Expected outputs:
- Downsampled BAM files (*.1000x.bam) in `./data/mapped_reads`
- Index stats files (*.1000x.idxstats) in `./data/mapped_reads`

4. Run SCiMS on downsampled simulated data:
```bash
bash 01_simulation/04_simulation_scims.sh
```

5. Run Rx, Ry, and BeXY on downsampled simulated data:
Execute the following commands to run Rx, Ry:
```bash
bash 01_simulation/05_simulation_rxry.sh
```

Execute the following commands to run BeXY:

BeXY installation guide for Linux: (original intructions can be found [Bexy](https://bitbucket.org/wegmannlab/bexy/wiki/Installation))
```bash
conda create -n bexy
conda activate bexy
conda install -n bexy -f ./envs/bexy.yaml

git clone https://bitbucket.org/WegmannLab/bexy.git
cd bexy
bash compile_bexy.sh

cp ./build/bexy "$CONDA_PREFIX/bin"
chmod +x "$CONDA_PREFIX/bin/bexy"
```

```bash
conda activate bexy

bash 01_simulation/06_simulation_bexy.sh
```

```bash
# generate bexy output in R
Rscript ./scripts/bexy.R ./01_simulation/bexy_output
```

### Part 2: Human Microbiome Project (HMP) dataset

1. Download the HMP dataset and map reads to the reference genomes:
```bash
snakemake -s 02_hmp/01_hmp_data_processing.smk --cores 4 --use-conda
```

2. Run SCiMS on HMP dataset:
```bash
snakemake -s 02_hmp/02_hmp_scims.smk --cores 4 --use-conda
```

3. Run Rx, Ry, and BeXY on HMP dataset:
Execute the following commands to run Rx, Ry:
```bash
bash 02_hmp/03_hmp_rxry.sh
``` 

Execute the following commands to run BeXY:
```bash
bash 02_hmp/04_hmp_bexy.sh
```

```bash
# generate bexy output in R
Rscript ./scripts/bexy.R ./02_hmp/results/bexy
```

### Part 3: Mouse Metagenomic Dataset

1. Download the mouse metagenomic dataset and map reads to the reference genomes:
```bash
snakemake -s 03_mouse/01_mouse_data_processing.smk --cores 4 --use-conda
```

2. Run SCiMS on mouse metagenomic dataset:
```bash
snakemake -s 03_mouse/02_mouse_scims.smk --cores 4 --use-conda
```

3. Run Rx, Ry, and BeXY on mouse metagenomic dataset:
Execute the following commands to run Rx, Ry:
```bash
bash 03_mouse/03_mouse_rxry.sh
``` 

Execute the following commands to run BeXY:
```bash
bash 03_mouse/04_mouse_bexy.sh
``` 

```bash
# generate bexy output in R
Rscript ./scripts/bexy.R ./03_mouse/results/bexy
```

### Part 4: Chicken Metagenomic Dataset

1. Download the chicken metagenomic dataset and map reads to the reference genomes:
```bash
snakemake -s 04_chicken/01_chicken_data_processing.smk --cores 4 --use-conda
```

2. Run SCiMS on chicken metagenomic dataset:
```bash
snakemake -s 04_chicken/02_chicken_scims.smk --cores 4 --use-conda
```

3. Run Rx, Ry, and BeXY on chicken metagenomic dataset:
Execute the following commands to run Rx, Ry:
```bash
bash 04_chicken/03_chicken_rxry.sh
```

Execute the following commands to run BeXY:
```bash
bash 04_chicken/04_chicken_bexy.sh
```

```bash
# generate bexy output in R
Rscript ./scripts/bexy.R
```

### Part 5: Generate figures
1. Generate figure 2:
```bash
python3 ./scripts/figure2.py
```

2. Generate figure 3:
```bash
python3 ./scripts/figure3.py
```

3. Generate figure 4:
```bash
python3 ./scripts/figure4.py

```

