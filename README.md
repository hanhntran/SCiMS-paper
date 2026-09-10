## SCiMS-paper

This repository contains reproducible scripts for the SCiMS manuscript.

This repository contains the pipelines used to generate
every analysis and figure in the paper. The tool itself lives at
[davenport-lab/SCiMS](https://github.com/davenport-lab/SCiMS)

> **Preprint / paper:** [https://www.biorxiv.org/content/10.64898/2026.02.17.705110v2.full](https://www.biorxiv.org/content/10.64898/2026.02.17.705110v2.full)

---

## Repository layout


| Directory                 | Contents                                                                           |
| ------------------------- | ---------------------------------------------------------------------------------- |
| `01_simulation/`          | CAMI-based simulations: read generation, mapping, downsampling, and all four tools |
| `02_hmp/`                 | Human Microbiome Project (four body sites), raw reads via dbGaP controlled access  |
| `02a_indian_metagenomes/` | Indian gut metagenomes (PRJNA397112), raw deposit                                  |
| `02b_hadza_metagenomes/`  | Hadza gut metagenomes (PRJEB49206), host-depleted deposit                          |
| `03_mouse/`               | Mouse metagenomes                                                                  |
| `04_chicken/`             | Chicken metagenomes (ZW sex-determination system)                                  |
| `05_baboon/`              | Baboon metagenomes                                                                 |
| `06_black_rhino/`         | Black rhino metagenomes                                                            |
| `07_cow/`                 | Cattle metagenomes                                                                 |
| `08_mesquite_lizard/`     | Mesquite lizard metagenomes                                                        |
| `09_pig/`                 | Pig metagenomes                                                                    |
| `data/`                   | Reference genomes, metadata tables, and intermediate outputs                       |
| `envs/`                   | Conda environment specifications                                                   |
| `figure_scripts/`         | Scripts that generate the main and supplementary figures                           |




### Environment setup

1. Clone the repository

```bash
git clone https://github.com/hanhntran/SCiMS-paper.git
cd SCiMS-paper
```

1. Install mamba: if you haven't installed mamba yet, use the following command to install it:

```bash
conda install -y -c conda-forge -c bioconda mamba
```

1. Create and activate the conda environment:

```bash
mamba env create -n scims-env -f ./envs/environment.yaml
mamba activate scims-env
```



#### BeXY

BeXY is compiled separately (original instructions:
[BeXY installation](https://bitbucket.org/wegmannlab/bexy/wiki/Installation)):

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



### Part 1: Simulation

1. Generate simulated reads:

```bash
bash 01_simulation/01_generate_host_reads.sh
```

Outputs: simulated male and female reference genomes in `./data/simulated_ref`,
simulated FASTQ reads in `./data/simulated_reads`.

1. Downsample host reads into different read depths:

```bash
bash 01_simulation/02_downsample_reads.sh
```

1. Mix in microbial reads from CAMI simulation oral sample:

```bash
bash 01_simulation/03_mix_host_cami.sh
```

1. Map host reads back to human genome

```bash
bash 01_simulation/04_map_mixed_reads.sh
```

1. Run the four tools:

```bash
bash 01_simulation/scims_hg38.sh     # SCiMS
bash 01_simulation/05_simulation_rxry.sh      # Rx and Ry

conda activate bexy
bash 01_simulation/run_bexy.sh      # BeXY
Rscript ./scripts/bexy.R ./01_simulation/bexy_output
```



## Part 2: Human cohorts



### 2a. Human Microbiome Project

Requires approved dbGaP access (phs000228); raw reads are used prior to host-read
screening.

```bash
bash 02_hmp/01_download_sra.sh
bash 02_hmp/02_map_reads.sh
bash 02_hmp/02_hmp_scims
bash 02_hmp/03_hmp_rxry.sh
bash 02_hmp/04_hmp_bexy.sh
Rscript ./scripts/bexy.R ./02_hmp/results/bexy
```



### 2b. Indian gut metagenomes (PRJNA397112)

```bash
bash 02a_indian_metagenomes/01_download_sra.sh
bash 02a_indian_metagenomes/02_map_reads.sh
bash 02a_indian_metagenomes/02_scims.sh
bash 02a_indian_metagenomes/03_rxry.sh
bash 02a_indian_metagenomes/04_bexy.sh
Rscript ./scripts/bexy.R ./02a_indian_metagenomes/results/bexy
```



### 2c. Hadza gut metagenomes (PRJEB49206)

The public deposit for this cohort is host-depleted (see
[note above](#a-note-on-host-read-filtering)). It is analyzed to characterize the
effect of host-read removal on sex inference, not as an intended-use benchmark.

```bash
bash 02b_hadza_metagenomes/01_download_sra.sh
bash 02b_hadza_metagenomes/02_map_reads.sh
bash 02b_hadza_metagenomes/02_scims.sh
bash 02b_hadza_metagenomes/03_rxry.sh
bash 02b_hadza_metagenomes/04_bexy.sh
Rscript ./scripts/bexy.R ./02b_hadza_metagenomes/results/bexy
```

---



## Part 3: Non-human hosts

The same four-step pattern applies to each species. Chicken uses a ZW
sex-determination system; species without a chromosome-level reference assembly with
annotated sex chromosomes were mapped to a closely related species' genome (see
Table 1 in the manuscript).

```bash
# replace <dir> and <name> with the dataset below
bash <dir>/01_download_sra.sh
bash <dir>/02_map_reads.sh
bash <dir>/02_scims.sh
bash <dir>/03_rxry.sh
bash <dir>/04_bexy.sh
Rscript ./scripts/bexy.R ./<dir>/results/bexy
```


| `<dir>`              | `<name>`          | Species                |
| -------------------- | ----------------- | ---------------------- |
| `03_mouse`           | `mouse`           | *Mus musculus*         |
| `04_chicken`         | `chicken`         | *Gallus gallus* (ZW)   |
| `05_baboon`          | `baboon`          | *Papio* spp.           |
| `06_black_rhino`     | `black_rhino`     | *Diceros bicornis*     |
| `07_cow`             | `cow`             | *Bos taurus*           |
| `08_mesquite_lizard` | `mesquite_lizard` | *Sceloporus grammicus* |
| `09_pig`             | `pig`             | *Sus domesticus*       |


---



## Part 4: Figures

```bash
python3 ./figure_scripts/figure2.py    # simulation benchmarking
python3 ./figure_scripts/figure3.py    # human cohort benchmarking
python3 ./figure_scripts/figure4.py    # cross-species benchmarking
```

Supplementary figures:

```bash
python3 ./figure_scripts/supp_fig_yfrac_by_sex_fecal.py \
    --hmp   data/dbGap_metadata_scims_updated_fecal.txt \
    --india data/indian_metadata_scims_updated.txt \
    --hadza data/hadza_PRJEB49206_metadata_scims_updated_filt.txt \
    --out   figures/Fig_S3_yfrac_by_sex_fecal
```

---

