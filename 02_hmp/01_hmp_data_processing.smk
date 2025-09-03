import os

main_dir = os.getcwd()

genome_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14"
genome_file = "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"


# Paths and files
REF_GENOME = os.path.join(main_dir, "data/ref_genome", "GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
BOWTIE2_INDEX_PREFIX = os.path.join(main_dir, "data/ref_genome", "GRCh38")
SCAFFOLD_FILE = os.path.join(main_dir, "data/ref_genome", "GRCh38_scaffold.txt")
SRA_FILE = os.path.join(main_dir, "data", "human_sra.txt")

raw_reads_dir = os.path.join(main_dir, "02_hmp/raw_reads")
trimmed_reads_dir = os.path.join(main_dir, "02_hmp/trimmed_reads")
mapped_reads_dir = os.path.join(main_dir, "02_hmp/mapped_reads")
ref_dir = os.path.join(main_dir, "data/ref_genome")

os.makedirs(raw_reads_dir, exist_ok=True)
os.makedirs(trimmed_reads_dir, exist_ok=True)
os.makedirs(mapped_reads_dir, exist_ok=True)
os.makedirs(ref_dir, exist_ok=True)

# Read SRA IDs
with open(SRA_FILE, 'r') as f:
    sra = [line.strip() for line in f]

rule all:
    input:
        expand(f"{mapped_reads_dir}/{{sra}}.idxstats", sra=sra),
        expand(BOWTIE2_INDEX_PREFIX + "{ext}", ext=[".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]),
        REF_GENOME

rule download_ref_genome:
    output:
        REF_GENOME
    shell:
        "wget {genome_url}/{genome_file} -O {output}"

rule index_ref_genome:
    input:
        REF_GENOME
    output:
        BOWTIE2_INDEX_PREFIX + ".1.bt2"
    conda:
        f"{main_dir}/envs/bowtie2_samtools.yaml"
    shell:
        "bowtie2-build {input} {BOWTIE2_INDEX_PREFIX}"

rule download_sra_reads:
    output:
        forward = f"{raw_reads_dir}/{{sra}}_1.fastq.gz",
        reverse = f"{raw_reads_dir}/{{sra}}_2.fastq.gz"
    params:
        sra_id = "{sra}"
    conda:
        f"{main_dir}/envs/sra-tools.yaml"
    shell:
        "fastq-dump --split-files --gzip --outdir {raw_reads_dir} {params.sra_id}"

rule trim_reads:
    input:
        forward = f"{raw_reads_dir}/{{sra}}_1.fastq.gz",
        reverse = f"{raw_reads_dir}/{{sra}}_2.fastq.gz"
    output:
        forward = f"{trimmed_reads_dir}/{{sra}}_1.trimmed.fastq.gz",
        reverse = f"{trimmed_reads_dir}/{{sra}}_2.trimmed.fastq.gz"
    conda:
        f"{main_dir}/envs/fastp.yaml"
    shell:
        """
        fastp -i {input.forward} -I {input.reverse} \
              -o {output.forward} -O {output.reverse} \
              -h {trimmed_reads_dir}/{wildcards.sra}.html \
              -j {trimmed_reads_dir}/{wildcards.sra}.json
        """

rule map_reads:
    input:
        forward = f"{trimmed_reads_dir}/{{sra}}_1.trimmed.fastq.gz",
        reverse = f"{trimmed_reads_dir}/{{sra}}_2.trimmed.fastq.gz",
        index = BOWTIE2_INDEX_PREFIX + ".1.bt2"
    output:
        bam = f"{mapped_reads_dir}/{{sra}}.sorted.bam"
    conda:
        f"{main_dir}/envs/bowtie2_samtools.yaml"
    shell:
        """
        bowtie2 -x {BOWTIE2_INDEX_PREFIX} \
                -1 {input.forward} -2 {input.reverse} | \
        samtools view -b -q 30 | \
        samtools sort -o {output.bam}
        """
        
rule add_read_groups:
    """
    Add read groups to the sorted BAM file. Read groups (RG) are required by Picard mark duplicates (next rule). 
    Because these reads are simulated, I just made up some random values to satisfy the RG requirements.
    """
    input:
        sorted_bam=f"{mapped_reads_dir}/{{sra}}.sorted.bam"
    output:
        rg_bam=f"{mapped_reads_dir}/{{sra}}.sorted.rg.bam"
    params:
        RGID="{sra}",
        RGLB="lib1",
        RGPL="illumina",
        RGPU="unit1",
        RGSM="{sra}"
    conda:
        f"{main_dir}/envs/picard.yaml"
    shell:
        """
        picard AddOrReplaceReadGroups \
            I={input.sorted_bam} \
            O={output.rg_bam} \
            RGID={params.RGID} \
            RGLB={params.RGLB} \
            RGPL={params.RGPL} \
            RGPU={params.RGPU} \
            RGSM={params.RGSM} \
            CREATE_INDEX=true
        """

rule remove_duplicates:
    input:
        bam = f"{mapped_reads_dir}/{{sra}}.sorted.rg.bam"
    output:
        bam = f"{mapped_reads_dir}/{{sra}}.rmdup.bam",
        metrics = f"{mapped_reads_dir}/{{sra}}.metrics.txt"
    conda:
        f"{main_dir}/envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            CREATE_INDEX=true
        """

rule idxstats:
    input:
        bam = f"{mapped_reads_dir}/{{sra}}.rmdup.bam"
    output:
        idxstats = f"{mapped_reads_dir}/{{sra}}.idxstats"
    conda:
        f"{main_dir}/envs/bowtie2_samtools.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats}
        """
