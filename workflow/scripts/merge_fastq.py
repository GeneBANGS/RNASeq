from subprocess import call

if len(snakemake.input) > 1:
    cmd = "cat"
    for i in snakemake.input:
        cmd=cmd+" "+i
    cmd=cmd+" > "
    cmd= cmd+snakemake.output[0]
    call(cmd, shell=True)
else:
    cmd = "cp "+snakemake.input[0]
    cmd=cmd+" "+snakemake.output[0]
    call(cmd, shell=True)