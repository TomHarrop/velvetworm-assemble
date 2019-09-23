#!/usr/bin/env python3

import multiprocessing
import pathlib2
import tempfile


#############
# FUNCTIONS #
#############

def busco_resolver(wildcards):
    return({
        'fasta': busco_inputs[wildcards.name]
        })


def resolve_path(x):
    return(str(pathlib2.Path(x).resolve()))


###########
# GLOBALS #
###########

ont_raw = 'data/all_passed_reads.fastq'

busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
flye_container = 'shub://TomHarrop/singularity-containers:flye_2.5'
minimap_container = 'shub://TomHarrop/singularity-containers:minimap2_2.11r797'
racon_container = 'library://tomharrop/default/genomics:racon_1.4.7'


########
# MAIN #
########

busco_inputs = {
    'flye': 'output/010_flye/assembly.fasta',
    'polish_lr': 'output/020_long_read_polishing/racon.fasta'
}

busco_lineages = [
    'arthropoda_odb9',
    'nematoda_odb9',
    'metazoa_odb9']


#########
# RULES #
#########

rule target:
    input:
        'output/020_long_read_polishing/aln.sam',
        expand(('output/099_busco/{lineage}/run_{assembly}/'
                'full_table_{assembly}.tsv'),
               lineage=busco_lineages,
               assembly=list(busco_inputs.keys()))


rule polish_long_reads:
    input:
        fasta = 'output/010_flye/assembly.fasta',
        aln = 'output/020_long_read_polishing/aln.sam',
        fq = ont_raw,
    output:
        'output/020_long_read_polishing/racon.fasta'
    log:
        'output/logs/020_long_read_polishing/racon.log'
    threads:
        multiprocessing.cpu_count()
    priority:
        0
    singularity:
        racon_container
    shell:
        'racon '
        '-t {threads} '
        '{input.fq} '
        '{input.aln} '
        '{input.fasta} '
        '> {output} '
        '2> {log}'


rule map_ont_reads:
    input:
        fasta = 'output/010_flye/assembly.fasta',
        fq = ont_raw
    output:
        'output/020_long_read_polishing/aln.sam'
    log:
        'output/logs/020_long_read_polishing/minimap.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        minimap_container
    shell:
        'minimap2 '
        '-t {threads} '
        '-ax '
        'map-ont '
        '{input.fasta} '
        '{input.fq} '
        '> {output} '
        '2> {log}'


rule flye:
    input:
        fq = ont_raw
    output:
        'output/010_flye/assembly.fasta'
    params:
        outdir = 'output/010_flye',
        size = '2g'
    threads:
        multiprocessing.cpu_count()
    priority:
        10
    log:
        'output/logs/010_flye.log'
    singularity:
        flye_container
    shell:
        'flye '
        '--iterations 1 '
        '--resume '
        '--nano-raw {input.fq} '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&>> {log}'

# generic BUSCO rule
rule busco:
    input:
        unpack(busco_resolver),
        lineage = 'data/busco/{lineage}'
    output:
        'output/099_busco/{lineage}/run_{name}/full_table_{name}.tsv'
    log:
        resolve_path('output/logs/099_busco/{lineage}-{name}.log')
    params:
        wd = 'output/099_busco/{lineage}',
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        lineage = lambda wildcards, input: resolve_path(input.lineage),
        tmpdir = tempfile.mkdtemp()
    threads:
        multiprocessing.cpu_count()
    priority:
        1
    singularity:
        busco_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'run_BUSCO.py '
        '--force '
        '--tmp_path {params.tmpdir} '
        '--in {params.fasta} '
        '--out {wildcards.name} '
        '--lineage {params.lineage} '
        '--cpu {threads} '
        '--mode genome '
        '&> {log}'

