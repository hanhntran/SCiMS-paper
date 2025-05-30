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

5. Run Rx, Ry, and BeXY on downsampled simulated data:
```bash
bash 01_simulation/05_simulation_rxry.sh
bash 01_simulation/06_simulation_bexy.sh
```

7. Generate the Figure 2:
```bash
python3 ./scripts/figure2.py
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
```bash
bash 02_hmp/03_hmp_rxry.sh
bash 02_hmp/04_hmp_bexy.sh
```

4. Generate the Figure 3:
```bash
python3 ./scripts/figure3_A.py
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
```bash
bash 03_mouse/03_mouse_rxry.sh
bash 03_mouse/04_mouse_bexy.sh
``` 

4. Generate the Figure 4:
```bash
python3 ./scripts/figure3_B.py
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
```bash
bash 04_chicken/03_chicken_rxry.sh
bash 04_chicken/04_chicken_bexy.sh
```

4. Generate the Figure 5:
```bash
python3 ./scripts/figure3_C.py  
```


