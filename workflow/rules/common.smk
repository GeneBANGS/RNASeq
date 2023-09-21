#######################
import errno
import pandas as pd
import os
import multiprocessing
from snakemake.utils import validate


report: "../report/workflow.rst"


# validate(config, schema="../schemas/config.schema.yaml")


samples = pd.read_table(config.get("samples"), index_col="sample")
units = pd.read_table(config.get("units"), index_col=["unit"], dtype=str)


## pipeline-related functions
def get_unit_fastqs(wildcards, samples, label="units", read_pair="fq"):
    for unit_set in samples.loc[wildcards.sample, [label]]:
        wildcards.sample
    return [
        expand_filepath(units.loc[x, [read_pair]].dropna()[0])
        for x in unit_set.split(",")
    ]


def get_odp(wildcards, samples, optical_dup="odp"):
    return "OPTICAL_DUPLICATE_PIXEL_DISTANCE={}".format(
        samples.loc[wildcards.sample, [optical_dup]].dropna()[0]
    )


## filepath functions
def resolve_results_filepath(basepath, outname):
    return os.path.join(basepath, outname)


def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT,
            os.strerror(errno.ENOENT) + " (path must be absolute)",
            filepath,
        )
    return filepath


def resolve_single_filepath(basepath, filename):
    return os.path.join(basepath, filename)


def resolve_multi_filepath(basepath, dictionary):
    for k, v in dictionary.items():
        dictionary[k] = os.path.join(basepath, v)
    return dictionary


## functions for system resources
def cpu_count():
    return multiprocessing.cpu_count()


def conservative_cpu_count(reserve_cores=1, max_cores=5):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)


def threads_calculator(read_type="pe", max_cores=99):
    if read_type == "se":
        if conservative_cpu_count(reserve_cores=1, max_cores=max_cores) > 2:
            return conservative_cpu_count(reserve_cores=1, max_cores=max_cores) / 2
        else:
            return 1
    else:
        if conservative_cpu_count(reserve_cores=1, max_cores=max_cores) > 4:
            return conservative_cpu_count(reserve_cores=1, max_cores=max_cores) / 4
        else:
            return 1
