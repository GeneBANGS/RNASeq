rule kallisto_index:
    input:
        resolve_single_filepath(
            config.get("resources").get("reference_path"),
            config.get("resources").get("transcriptome_fasta"),
        ),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "kallisto/index/transcriptome.kidx.transcripts",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/kallisto.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=1, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Build Transcriptome Index for Kallisto with {threads} threads on the following fasta file: {input}."
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "logs/kallisto/index_kallisto.log"
        ),
    shell:
        "kallisto index "
        "-i {output} "
        "{input} "
        ">& {log} "


rule kallisto:
    input:
        r1=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/trimmed/{sample}-R1-trimmed.fq.gz",
        ),
        r2=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/trimmed/{sample}-R2-trimmed.fq.gz",
        ),
        index=rules.kallisto_index.output,
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "kallisto/{sample}/abundance.h5"
        ),
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "kallisto/{sample}/abundance.tsv"
        ),
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "kallisto/{sample}/run_info.json"
        ),
    params:
        outdir=lambda w, output: os.path.split(output[0])[0],
        params=config.get("params").get("kallisto").get("arguments"),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/kallisto.yaml"
        )
    threads:
        conservative_cpu_count(
            reserve_cores=1, max_cores=int(config.get("resources").get("max_cores"))
        )
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Kallisto Transcript Quantification with {threads} threads on the following fastq files: {input.r1}/{input.r2}."
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/kallisto/{sample}.kallisto.log",
        ),
    shell:
        "kallisto quant "
        "-i {input.index} "
        "-o {params.outdir} "
        "{params.params} "
        "-t {threads} "
        "{input.r1} {input.r2} "
        ">& {log}"
