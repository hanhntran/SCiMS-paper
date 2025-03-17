import os

# Directories and paths
dir_path = os.getcwd()
ref_genome_index = os.path.join(dir_path, "data/ref_genome/GRCh38_latest_genomic.fna")
simulated_reads_dir = os.path.join(dir_path, "data/simulated_reads")
map_dir = os.path.join(dir_path, "data/mapped_reads")

# Create output directory
os.makedirs(map_dir, exist_ok=True)

rule all:
    input:
        expand(f"{map_dir}/S{{sex}}100000000.sorted.rmdup.bam", sex=["M", "F"])

rule map_simulated_reads:
    input:
        forward_reads=f"{dir_path}/data/simulated_reads/S{{sex}}100000000_1.fastq.gz",
        reverse_reads=f"{dir_path}/data/simulated_reads/S{{sex}}100000000_2.fastq.gz",
        index=ref_genome_index + ".bwt"
    output:
        sorted_bam=f"{dir_path}/data/mapped_reads/S{{sex}}100000000.sorted.bam"
    params:
        index=ref_genome_index
    conda:
        "./envs/bowtie2_samtools.yaml"
    shell:
        """
        bowtie2 -x {params.index} -1 {input.forward_reads} -2 {input.reverse_reads} | \
        samtools view -b -q 30 | \
        samtools sort -o {output}
        """

rule add_read_groups:
    input:
        sorted_bam=f"{map_dir}/S{{sex}}100000000.sorted.bam"
    output:
        rg_bam=f"{map_dir}/S{{sex}}100000000.sorted.rg.bam"
    params:
        RGID="S{sex}100000000",
        RGLB="lib1",
        RGPL="illumina",
        RGPU="unit1",
        RGSM="S{sex}100000000"
    conda:
        "./envs/picard.yaml"
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

rule mark_duplicates:
    input:
        rg_bam=f"{map_dir}/S{{sex}}100000000.sorted.rg.bam"
    output:
        rmdup_bam=f"{map_dir}/S{{sex}}100000000.sorted.rmdup.bam",
        metrics=f"{map_dir}/S{{sex}}100000000.duplication_metrics.txt"
    conda:
        "./envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates \
            I={input.rg_bam} \
            O={output.rg_bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            CREATE_INDEX=true
        """

