rule star_index:
    input:
        resolve_single_filepath(
            config.get("resources").get("reference_path"),
            config.get("resources").get("reference_fa"),
        ),
    output:
        length=ensure(
            resolve_results_filepath(
                config.get("paths").get("results_dir"), "star/index/chrLength.txt"
            ),
            non_empty=True,
        ),
    params:
        gtf=resolve_single_filepath(
            config.get("resources").get("reference_path"),
            config.get("resources").get("gtf"),
        ),
        genomeDir=lambda w, output: os.path.split(output[0])[0],
        overhang=config.get("params").get("star").get("overhang"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "logs/star/star_index.log"
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/star.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Indexing Genome with {threads} threads on the following fasta file: {input}."
    shell:
        "STAR "
        "--runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeFastaFiles {input} "
        "--sjdbGTFfile {params.gtf} "
        "--sjdbOverhang {params.overhang} "
        "--outTmpDir {resources.tmpdir}/STARtmp "
        "--genomeDir {params.genomeDir} "
        "&> {log} "


rule star_2pass_mapping:
    input:
        read1=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/trimmed/{sample}-R1-trimmed.fq.gz",
        ),
        read2=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/trimmed/{sample}-R2-trimmed.fq.gz",
        ),
        index=rules.star_index.output.length,
    output:
        bam=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
            )
        ),
        log=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "star/{sample}/{sample}.Log.final.out",
        ),
    params:
        genomedir=lambda w, input: os.path.split(input.index)[0],
        sample="{sample}",
        platform=config.get("params").get("star").get("platform"),
        center=config.get("params").get("star").get("sequencing_center"),
        out_basename=lambda w, output: output[0].split(".", 1),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/star/{sample}/{sample}_star_map.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/star.yaml"
        )
    threads:
        conservative_cpu_count(
            reserve_cores=1, max_cores=int(config.get("resources").get("max_cores"))
        )
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=30000,
    message:
        "STAR RNA-Seq Aligner 2-pass mapping with {threads} threads and {resources.mem_mb} Mbytes for the following file: {input.read1} and {input.read2}."
    shell:
        "STAR "
        "--runMode alignReads "
        "--genomeDir {params.genomedir} "
        "--runThreadN {threads} "
        "--readFilesIn {input.read1} {input.read2} "
        " --outSAMattrRGline  ID:{params.sample} SM:{params.sample} LB:library PL:{params.platform}  PU:{params.platform} CN:{params.center} "
        "--readFilesCommand zcat "
        "--twopass1readsN -1 "
        "--outStd Log "
        "--outSAMtype BAM SortedByCoordinate "
        "--twopassMode Basic "
        "--limitBAMsortRAM 30000000000 "
        "--outSAMmapqUnique 60 "
        "--outTmpDir {resources.tmpdir}/STARtmp "
        "--outFileNamePrefix {params.out_basename} "
        "&> {log} "


rule move_bam_files:
    input:
        bam=rules.star_2pass_mapping.output.bam,
    output:
        bam=protected(
            resolve_results_filepath(
                config.get("paths").get("results_dir"), "reads/bam/{sample}.bam"
            )
        ),
    params:
        sample="{sample}",
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/star/{sample}/{sample}_star_map_move.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/star.yaml"
        )
    threads:
        conservative_cpu_count(
            reserve_cores=1, max_cores=int(config.get("resources").get("max_cores"))
        )
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
        mem_mb=30000,
    message:
        "Rename and move with {threads} threads and {resources.mem_mb} Mbytes for the following BAM file: {input.bam}."
    shell:
        "mv "
        "{input.bam} {output.bam} "
        "&> {log} "


rule index_bam:
    input:
        bam=rules.move_bam_files.output.bam,
    output:
        bai=protected(
            resolve_results_filepath(
                config.get("paths").get("results_dir"), "reads/bam/{sample}.bam.bai"
            )
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    threads:
        conservative_cpu_count(
            reserve_cores=1, max_cores=int(config.get("resources").get("max_cores"))
        )
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    message:
        "Build Index for STAR mapped reads with {threads} threads on the following BAM file: {input.bam}."
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "logs/samtools/{sample}.index.log"
        ),
    shell:
        "samtools index "
        "{input.bam} "
        ">& {log} "
