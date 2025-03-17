import os

dir_path = os.getcwd()

# Define paths and files
genome_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14"
genome_file = "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"

ref_dir = os.path.join(dir_path, "data/ref_genome")
simulated_ref_dir = os.path.join(dir_path, "data/simulated_ref")
simulated_reads_dir = os.path.join(dir_path, "data/simulated_reads")

os.makedirs(ref_dir, exist_ok=True)
os.makedirs(simulated_ref_dir, exist_ok=True)
os.makedirs(simulated_reads_dir, exist_ok=True)

rule all:
    input:
        expand(f"{simulated_reads_dir}/S{{sex_abbr}}100000000_{{pair}}.fastq.gz", 
               sex_abbr=["M","F"], pair=[1,2])

rule download_ref_genomes:
    output:
        f"{ref_dir}/{genome_file}"
    conda: "./envs/python3.9.yaml"
    params:
        genome_url=genome_url,
        genome_file=genome_file
    shell:
        "wget {params.genome_url}/{params.genome_file} -O {output}"

rule create_ref_genomes_for_male_and_female:
    input:
        ref=f"{ref_dir}/{genome_file}"
    output:
        male=f"{simulated_ref_dir}/simulated_ref_male.fasta",
        female=f"{simulated_ref_dir}/simulated_ref_female.fasta"
    conda: "./envs/python3.9.yaml"
    shell:
        "python ./scripts/create_ref.py -r {input.ref} -o {simulated_ref_dir}/simulated_ref"

rule simulate_reads_male:
    input: rules.create_ref_genomes_for_male_and_female.output.male
    output:
        fastq1=f"{simulated_reads_dir}/SM100000000_1.fastq.gz",
        fastq2=f"{simulated_reads_dir}/SM100000000_2.fastq.gz"
    conda: "./envs/wgsim-1.0.yaml"
    shell: 
        """
        wgsim -N 100000000 -1 150 -2 150 {input} male_1.fastq male_2.fastq
        gzip male_1.fastq male_2.fastq
        mv male_1.fastq.gz {output.fastq1}
        mv male_2.fastq.gz {output.fastq2}
        """

rule simulate_reads_female:
    input: rules.create_ref_genomes_for_male_and_female.output.female
    output:
        fastq1=f"{simulated_reads_dir}/SF100000000_1.fastq.gz",
        fastq2=f"{simulated_reads_dir}/SF100000000_2.fastq.gz"
    conda: "./envs/wgsim-1.0.yaml"
    shell:  
        """
        wgsim -N 100000000 -1 150 -2 150 {input} female_1.fastq female_2.fastq
        gzip female_1.fastq female_2.fastq
        mv female_1.fastq.gz {output.fastq1}
        mv female_2.fastq.gz {output.fastq2}
        """ 
