
rule fastqc:
    input:
        rules.merge_read1_fastq.output,
        rules.merge_read2_fastq.output
    output:
        html_r1 = resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/{sample}-R1.html"),
        zip_r1 = resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/{sample}-R1_fastqc.zip"),
        html_r2= resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/{sample}-R2.html"),
        zip_r2= resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc/{sample}-R2_fastqc.zip")
    params:
        outdir=lambda w, output: os.path.split(output[0])[0]
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/fastqc/{sample}.log")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/fastqc.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=1, max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=6000
    message: "."
    shell:
        "fastqc "
        "{input} "
        "--outdir {params.outdir} "
        "--quiet "
        ">& {log}"


rule fastqc_trimmed:
    input:
        rules.rename_trimmed_fastqs.output.read1,
        rules.rename_trimmed_fastqs.output.read2
    output:
        html_r1 = resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc_trimmed/{sample}-R1.html"),
        zip_r1 = resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc_trimmed/{sample}-R1_fastqc.zip"),
        html_r2= resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc_trimmed/{sample}-R2.html"),
        zip_r2= resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/fastqc_trimmed/{sample}-R2_fastqc.zip")
    params:
        outdir=lambda w, output: os.path.split(output[0])[0]
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/fastqc_trimmed/{sample}.log")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/fastqc.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=1, max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=6000
    message: "."
    shell:
        "fastqc "
        "{input} "
        "--outdir {params.outdir} "
        "--quiet "
        ">& {log}"


rule multiqc:
    input:
        fastqc_r1=expand(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "qc/fastqc/{sample.sample}-R1_fastqc.zip",
            ),
            sample=samples.reset_index().itertuples(),
        ),
        fastqc_r2=expand(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "qc/fastqc/{sample.sample}-R2_fastqc.zip",
            ),
            sample=samples.reset_index().itertuples(),
        ),
        fastqc_trimmed_r1=expand(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "qc/fastqc_trimmed/{sample.sample}-R1_fastqc.zip",
            ),
            sample=samples.reset_index().itertuples(),
        ),
        fastqc_trimmed_r2=expand(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "qc/fastqc_trimmed/{sample.sample}-R2_fastqc.zip",
            ),
            sample=samples.reset_index().itertuples(),
        ),
        rseqc=expand(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "qc/rseqc/{sample.sample}.rseqc_complete",
            ),
            sample=samples.reset_index().itertuples(),
        ),
        star=expand(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "star/{sample.sample}/{sample.sample}.Log.final.out",
            ),
            sample=samples.reset_index().itertuples(),
        ),
        kallisto=expand(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "logs/kallisto/{sample.sample}.kallisto.log",
            ),
            sample=samples.reset_index().itertuples(),
        ),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/multiqc.html")
    params:
        params=config.get("params").get("multiqc").get("arguments"),
        outdir="qc",
        outname="multiqc.html",
        fastqc=lambda w, input: os.path.split(input.fastqc_r1[0])[0],
        trimming=lambda w, input: os.path.split(input.fastqc_trimmed_r1[0])[0],
        rseqc=lambda w, input: os.path.split(input.rseqc[0])[0],
        star=lambda w, input: os.path.split(input.star[0])[0],
        kallisto=lambda w, input: os.path.split(input.kallisto[0])[0],
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/multiqc.yaml")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/multiqc/multiqc.log")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=6000
    message: "Generate MultiQC Report with {threads} threads and {resources.mem_mb} Mbytes."
    shell:
        "multiqc "
        "{params.fastqc} "
        "{params.trimming} "
        "{params.rseqc} "
        "{params.star} "
        "{params.kallisto} "
        "{params.params} "
        "-o {params.outdir} "
        "-n {params.outname} "
        ">& {log}"




