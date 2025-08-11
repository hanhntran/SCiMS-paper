import os

# Directories and paths
main_dir = os.getcwd()

genome_file = "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"

ref_dir = os.path.join(main_dir, "data/ref_genome")
ref_genome_index = os.path.join(ref_dir, genome_file)
simulated_reads_dir = os.path.join(main_dir, "01_simulation/simulated_reads")
map_dir = os.path.join(main_dir, "01_simulation/mapped_reads")

# Create output directory
os.makedirs(map_dir, exist_ok=True)

rule all:
    input:
        expand(f"{map_dir}/S{{sex}}100000000.sorted.rmdup.bam", sex=["M", "F"])

rule map_simulated_reads:
    """
    Map simulated reads to the reference genome using Bowtie2.
    """
    input:
        forward_reads=f"{simulated_reads_dir}/S{{sex}}100000000_1.fastq.gz",
        reverse_reads=f"{simulated_reads_dir}/S{{sex}}100000000_2.fastq.gz",
        index=ref_genome_index + ".bwt"
    output:
        sorted_bam=f"{map_dir}/S{{sex}}100000000.sorted.bam"
    params:
        index=ref_genome_index
    conda:
        f"{main_dir}/envs/bowtie2_samtools.yaml"
    shell:
        """
        bowtie2 -x {params.index} -1 {input.forward_reads} -2 {input.reverse_reads} | \
        samtools view -b -q 30 | \
        samtools sort -o {output}
        """

rule add_read_groups:
    """
    Add read groups to the sorted BAM file. Read groups (RG) are required by Picard mark duplicates (next rule). 
    Because these reads are simulated, I just made up some random values to satisfy the RG requirements.
    """
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

rule mark_duplicates:
    """
    Mark duplicates in the sorted BAM file
    """
    input:
        rg_bam=f"{map_dir}/S{{sex}}100000000.sorted.rg.bam"
    output:
        rmdup_bam=f"{map_dir}/S{{sex}}100000000.sorted.rmdup.bam",
        metrics=f"{map_dir}/S{{sex}}100000000.duplication_metrics.txt"
    conda:
        f"{main_dir}/envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates \
            I={input.rg_bam} \
            O={output.rmdup_bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            CREATE_INDEX=true
        """

