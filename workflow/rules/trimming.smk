rule merge_read1_fastq:
    input:
        lambda wildcards: get_unit_fastqs(wildcards, samples, read_pair='fq1')
    output:
        temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/untrimmed/merged/{sample}-R1.fq.gz"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/python_script.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3500
    script:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/scripts/merge_fastq.py"
        )


rule merge_read2_fastq:
    input:
        lambda wildcards: get_unit_fastqs(wildcards, samples, read_pair='fq2')
    output:
        temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/untrimmed/merged/{sample}-R2.fq.gz"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/python_script.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3500
    script:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/scripts/merge_fastq.py"
        )


rule trimming:
    input:
        read1=rules.merge_read1_fastq.output,
        read2=rules.merge_read2_fastq.output
    output:
        read1=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample}-R1_val_1.fq.gz")),
        read1_trimming_report=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample}-R1.fq.gz_trimming_report.txt"),
        read2=temp(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample}-R2_val_2.fq.gz")),
        read2_trimming_report=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample}-R2.fq.gz_trimming_report.txt")
    params:
        extra=config.get("params").get("trim_galore").get("arguments"),
        outdir=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/"),
        qc_dir=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/trim_galore/{sample}.log")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/trim_galore.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=0, max_cores=2)
    resources:
        tmpdir=config.get("paths").get("tmp_dir")
    message: "Trimming reads with TRIM_GALORE with {threads} threads for the following files {input.read1}{input.read2}."
    shell:
        "mkdir -p {params.qc_dir}; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"


rule rename_trimmed_fastqs:
    input:
        read1=rules.trimming.output.read1,
        read2=rules.trimming.output.read2
    output:
        read1=protected(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample}-R1-trimmed.fq.gz")),
        read2=protected(resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample}-R2-trimmed.fq.gz"))
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/bash/{sample}.log")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bash.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=1, max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir")
    message: "Rename fastq files {input.read1}{input.read2}."
    shell:
        "mv {input.read1} {output.read1} && "
        "mv {input.read2} {output.read2} "
        ">& {log}"