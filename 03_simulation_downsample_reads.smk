import os

# Directories and Parameters
dir_path = os.getcwd()
map_dir = os.path.join(dir_path, "data/mapped_reads")
depth = ["150", "200", "250", "300", "350", "400", "450", "1000", "5000", "10000", "100000", "1000000"]
sex = ["F", "M"]
iteration = list(range(1, 1001))

rule all:
    input:
        expand(os.path.join(map_dir, "S{iteration}{sex}{depth}.1000x.idxstats"), 
               iteration=iteration, sex=sex, depth=depth)

rule downsample_bam:
    input:
        bam=os.path.join(map_dir, "S1{sex}100000000.sorted.rmdup.bam")
    output:
        bam=os.path.join(map_dir, "S{iteration}{sex}{depth}.simulated.downsampled.1000x.bam")
    conda:
        "./envs/python_downsample.yaml"
    shell:
        """
        python3 scripts/finite_downsampler.py -t 157093441 -n {wildcards.depth} \
        -i {input.bam} -o {output.bam}
        """

rule index_downsampled_bam:
    input:
        bam=os.path.join(map_dir, "S{iteration}{sex}{depth}.simulated.downsampled.1000x.bam")
    output:
        bai=os.path.join(map_dir, "S{iteration}{sex}{depth}.simulated.downsampled.1000x.bam.bai")
    conda:
        "./envs/samtools1.19.yaml"
    shell:
        """
        samtools index {input.bam}
        """

rule generate_idxstats:
    input:
        bam=os.path.join(map_dir, "S{iteration}{sex}{depth}.simulated.downsampled.1000x.bam"),
        bai=os.path.join(map_dir, "S{iteration}{sex}{depth}.simulated.downsampled.1000x.bam.bai")
    output:
        idxstats=os.path.join(map_dir, "S{iteration}{sex}{depth}.1000x.idxstats")
    conda:
        "./envs/samtools1.19.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats}
        """

