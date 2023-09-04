rule bam_stat:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.bam_stat.txt")
    params:
        min_map_qual=config.get("params").get("rseqc").get("bamstat").get("min_map_qual")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=1, max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "Executing RSeQC BAMSTAT with {threads} threads on the following files {input.bam}."
    shell:
        "bam_stat.py "
        "-i {input.bam} "
        "-q {params.min_map_qual}"
        ">& {output}"


rule smatools_flagstat:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/samtools/{sample}.flagstat.txt")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=1, max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "Executing SAMTOOLS FLAGSTAT with {threads} threads on the following files {input.bam}."
    shell:
        "samtools flagstat "
        "-i {input.bam} "
        ">& {output}"


rule rseqc_read_distribution:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.read_distribution.txt")
    params:
#        ref=resolve_single_filepath(*references_abs_path(ref="rseqc_reference"), config.get("refseq")),
        ref=resolve_single_filepath(config.get("resources").get("rseqc_data"), config.get("resources").get("refseq_bed")),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=1, max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "Executing RSeQC READ DISTRIBUTION with {threads} threads on the following files {input.bam}."
    shell:
        "read_distribution.py "
        "-r {params.ref} "
        "-i {input.bam} "
        ">& {output}"


rule genebody_coverage:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.geneBodyCoverage.txt")
    params:
        out_basename=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}"),
        ref=resolve_single_filepath(config.get("resources").get("rseqc_data"), config.get("resources").get("housekeeping_genes"))
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/rseqc/genebody_coverage/{sample}_genebodycoverage.log")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1, max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "Executing RSeQC Gene Body Coverage with {threads} threads on the following files {input.bam}."
    shell:
        "geneBody_coverage.py "
        "-r {params.ref} "
        "-i {input.bam} "
        "-o {params.out_basename} "
        ">& {log} "

rule rseqc_junction_annotation:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.junction.txt")
    params:
        out_basename=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}"),
        ref=resolve_single_filepath(config.get("resources").get("rseqc_data"), config.get("resources").get("refseq_bed"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "."
    shell:
        "junction_annotation.py "
        "-r {params.ref} "
        "-i {input.bam} "
        "-o {params.out_basename} "
        "2> {output.out}"

rule rseqc_junction_saturation:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        plotr=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.junctionSaturation_plot.r")
    params:
        out_basename=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}"),
        ref=resolve_single_filepath(config.get("resources").get("rseqc_data"), config.get("resources").get("refseq_bed")),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.junctionSaturation.txt")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "."
    shell:
        "junction_saturation.py "
        "-r {params.ref} "
        "-i {input.bam} "
        "-o {params.out_basename} "
        "&> {log}"


rule rseqc_GC:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.GC.xls")
    params:
        out_basename=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "."
    shell:
        "read_GC.py "
        "-i {input.bam} "
        "-o {params.out_basename}"


rule rseqc_infer_experiment:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.infer_experiment.txt")
    params:
        ref=resolve_single_filepath(config.get("resources").get("rseqc_data"), config.get("resources").get("refseq_bed"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "."
    shell:
        "infer_experiment.py "
        "-r {params.ref} "
        "-i {input.bam} "
        ">& {output}"


rule rseqc_inner_distance:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.inner_distance.txt")
    params:
        out_basename = resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}"),
        ref=resolve_single_filepath(config.get("resources").get("rseqc_data"), config.get("resources").get("refseq_bed")),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "."
    shell:
        "inner_distance.py "
        "-r {params.ref} "
        "-i {input.bam} "
        "-o {params.out_basename}"


rule rseqc_read_duplication:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        out1=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.dup.pos.DupRate.xls")
    params:
        out_basename=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "."
    shell:
        "read_duplication.py "
        "-i {input.bam} "
        "-o {params.out_basename}"


rule rseqc_RPKM_saturation:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.saturation.pdf")
    params:
        out_basename=resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}"),
        ref=resolve_single_filepath(config.get("resources").get("rseqc_data"), config.get("resources").get("refseq_bed")),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/rseqc/{sample}.RPKM_saturation.log")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "."
    shell:
        "RPKM_saturation.py "
        "-r {params.ref} "
        "-i {input.bam} "
        "-o {params.out_basename} "
        ">& {log} "


rule rseqc_tin:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam"),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/bam/{sample}.bam.bai")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.tin.summary.txt")
    params:
        ref=resolve_single_filepath(config.get("resources").get("rseqc_data"), config.get("resources").get("refseq_bed"))
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "."
    shell:
        "tin.py "
        "-r {params.ref} "
        "-i {input.bam} "
        ">& {output} "


rule check_rseqc:
    input:
        rules.rseqc_tin.output,
        rules.rseqc_GC.output,
        rules.bam_stat.output,
        rules.rseqc_RPKM_saturation.output.out,
        rules.rseqc_inner_distance.output,
        rules.genebody_coverage.output,
        rules.rseqc_infer_experiment.output,
        rules.rseqc_junction_annotation.output.out,
        rules.rseqc_junction_saturation.output.plotr,
        rules.rseqc_read_distribution.output,
        rules.rseqc_read_duplication.output.out1
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"qc/rseqc/{sample}.rseqc_complete")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/rseqc.yaml")
    threads: conservative_cpu_count(reserve_cores=1,max_cores=int(config.get("resources").get("max_cores")))
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=3000
    message: "."
    shell:
        "touch {output} "