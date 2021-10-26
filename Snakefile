#!/usr/bin/env python3

import multiprocessing
import pathlib
import tempfile


#############
# FUNCTIONS #
#############

def busco_resolver(wildcards):
    return({
        'fasta': busco_inputs[wildcards.name]
        })


def resolve_path(x):
    return(str(pathlib.Path(x).resolve()))


###########
# GLOBALS #
###########

ont_raw = 'data/all_passed_reads.fastq'

bbmap_container = 'docker://ghcr.io/deardenlab/container-bbmap:bbmap_38.90'
biopython = 'shub://TomHarrop/singularity-containers:biopython_1.73'
busco = 'docker://ezlabgva/busco:v4.0.5_cv1'
flye_container = 'docker://ghcr.io/tomharrop/container-flye:2.9'
minimap_container = 'shub://TomHarrop/singularity-containers:minimap2_2.17r941'
pigz = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
porechop = 'shub://TomHarrop/ont-containers:porechop_0.2.4'
racon_chunks = 'shub://TomHarrop/racon-chunks:racon-chunks_0.0.5'
racon_container = ('docker://quay.io/tomharrop/genomics:'
                   'racon_ededb83-nvidia_410-bionic')
samtools_container = 'docker://ghcr.io/deardenlab/container-samtools:samtools_1.12'


########
# MAIN #
########

busco_inputs = {
    'flye': 'output/020_flye-polish/polished_1.fasta',
    'racon_lr': 'output/020_long_read_polishing/racon_lr.fasta',
    'racon_sr': 'output/030_short_read_polishing/racon_sr.fasta'
}

busco_lineages = [
    'arthropoda_odb9',
    # 'nematoda_odb9',
    'metazoa_odb9']

# chunkiness
n_chunks = 100
all_chunks = [str(x) for x in range(0, n_chunks)]
# all_chunks = [756, 817]


#########
# RULES #
#########

rule target:
    input:
        'output/030_short_read_polishing/racon_sr.fasta',
        expand(('output/099_busco/'
                '{name}/run_metazoa_odb10/'
                'full_table.tsv'),
               name=list(busco_inputs.keys()))


# short read racon chunks (no GPU)
rule racon_chunks:
    input:
        assembly = 'output/020_long_read_polishing/racon_lr.fasta',
        reads = 'output/030_short_read_polishing/short_reads.fq'
    output:
        'output/030_short_read_polishing/racon_sr.fasta'
    params:
        outdir = 'output/030_short_read_polishing',
        output_filename = 'racon_sr.fasta',
        chunks = '1000'
    log:
        'output/logs/030_short_read_polishing/racon_chunks.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        racon_chunks
    shell:
        'racon_chunks '
        '--reads {input.reads} '
        '--assembly {input.assembly} '
        '--outdir {params.outdir} '
        '--output_filename {params.output_filename} '
        '--threads {threads} '
        '--chunks {params.chunks} '
        '--wait_min 120 '
        '&> {log}'

rule combine_illumina:
    input:
        expand('output/030_short_read_polishing/run{run}.fastq',
               run=['1', '2'])
    output:
        'output/030_short_read_polishing/short_reads.fq'
    singularity:
        bbmap_container
    shell:
        'cat {input} > {output}'

rule trim_illumina:
    input:
        'output/030_short_read_polishing/run{run}_filter.fastq'
    output:
        p = pipe('output/030_short_read_polishing/run{run}.fastq'),
        stats = 'output/030_short_read_polishing/run{run}_trim.txt'
    log:
        'output/logs/030_short_read_polishing/run{run}_trim.log'
    params:
        trim = '/adapters.fa'
    threads:
        10
    singularity:
        bbmap_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.stats} '
        '>> {output.p} '
        '2> {log} '


rule filter_illumina:
    input:
        'output/030_short_read_polishing/repair_run{run}.fastq'
    output:
        p = pipe('output/030_short_read_polishing/run{run}_filter.fastq'),
        stats = 'output/030_short_read_polishing/run{run}_filter.txt'
    log:
        'output/logs/030_short_read_polishing/run{run}_filter.log'
    params:
        filter = '/phix174_ill.ref.fa.gz'
    threads:
        10
    singularity:
        bbmap_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.stats} '
        '>> {output.p} '
        '2> {log}'

rule repair_illumina:
    input:
        r1 = 'data/illumina_run{run}/CC481_R1.fq.gz',
        r2 = 'data/illumina_run{run}/CC481_R2.fq.gz',
    output:
        p = pipe('output/030_short_read_polishing/repair_run{run}.fastq'),
    log:
        'output/logs/030_short_read_polishing/repair_run{run}.log'
    singularity:
        bbmap_container
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        '>> {output.p} '
        '2> {log}'

# long read racon chunks
rule combine_racon_chunks:
    input:
        expand(('output/020_long_read_polishing/racon-chunks/'
                '{chunk}.fasta'),
               chunk=all_chunks)
    output:
        'output/020_long_read_polishing/racon_lr.fasta'
    log:
        'output/logs/020_long_read_polishing/combine_racon_chunks.log'
    singularity:
        racon_container
    shell:
        'cat {input} > {output} 2> {log}'    


rule polish_long_reads:
    input:
        fasta = 'output/015_genome-chunks/chunk_{chunk}.fasta',
        aln = 'output/020_long_read_polishing/bam-chunks/chunk_{chunk}.sam',
        fq = 'output/020_long_read_polishing/read-chunks/chunk_{chunk}.fq'
    output:
        temp('output/020_long_read_polishing/racon-chunks/{chunk}.fasta')
    params:
        wait_time = '60m',
        cuda_batches = 65
    log:
        'output/logs/020_long_read_polishing/racon-{chunk}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        racon_container
    shell:
        'timeout {params.wait_time} '
        'racon '
        '-t {threads} '
        '--cudapoa-batches {params.cuda_batches} '
        '{input.fq} '
        '{input.aln} '
        '{input.fasta} '
        '> {output} '
        '2> {log}'

# map long reads, chunk and polish
rule retrieve_reads:
    input:
        read_ids = expand(('output/020_long_read_polishing/read-ids/'
                           'chunk_{chunk}.txt'),
                          chunk=all_chunks),
        fastq = ont_raw
    output:
        temp(expand(('output/020_long_read_polishing/read-chunks/'
                     'chunk_{chunk}.fq'),
                    chunk=all_chunks))
    params:
        outdir = 'output/020_long_read_polishing/read-chunks/'
    log:
        'output/logs/020_long_read_polishing/retrieve_reads.log'
    singularity:
        biopython
    script:
        'src/retrieve_reads.py'

rule extract_read_ids:
    input:
        'output/020_long_read_polishing/bam-chunks/chunk_{chunk}.sam'
    output:
        'output/020_long_read_polishing/read-ids/chunk_{chunk}.txt'
    log:
        'output/logs/020_long_read_polishing/extract_read_ids_{chunk}.log'
    threads:
        1
    singularity:
        samtools_container
    shell:
        'samtools view  {input} '
        '| cut -f1 '
        '| sort '
        '| uniq '
        '> {output} '
        '2> {log}'

rule chunk_bam:
    input:
        bam = 'output/020_long_read_polishing/aln_sorted.bam',
        bai = 'output/020_long_read_polishing/aln_sorted.bam.bai',
        contig_list = 'output/015_genome-chunks/chunk_{chunk}_contigs.txt'
    output:
        temp('output/020_long_read_polishing/bam-chunks/chunk_{chunk}.sam')
    log:
        'output/logs/020_long_read_polishing/view_{chunk}.log'
    threads:
        1
    singularity:
        samtools_container
    shell:
        # the horrible sed command replaces the newlines with spaces in
        # contig_list. This allows the workflow to run before contig_list is
        # created. Adapted from https://stackoverflow.com/a/1252191
        'contigs="$(sed -e \':a\' -e \'N\' -e \'$!ba\' '
        '-e \'s/\\n/ /g\' {input.contig_list})" ; '
        'samtools view '
        '-h '
        '-F 2308 '  # read unmapped (0x4); not primary alignment (0x100);     supplementary alignment (0x800)
        '-O SAM '
        '{input.bam} '
        '${{contigs}} '
        '> {output} '
        '2> {log}'

rule index_bam:
    input:
        bam = 'output/020_long_read_polishing/aln_sorted.bam'
    output:
        bai = 'output/020_long_read_polishing/aln_sorted.bam.bai'
    log:
        'output/logs/020_long_read_polishing/index_bam.log'
    threads:
        1
    singularity:
        samtools_container
    shell:
        'samtools index {input.bam} {output.bai} '
        '2> {log}'

rule sort_sam:
    input:
        sam = 'output/020_long_read_polishing/aln.sam',
        fasta = 'output/010_flye/assembly.fasta'
    output:
        bam = 'output/020_long_read_polishing/aln_sorted.bam'
    log:
        'output/logs/020_long_read_polishing/sort_sam.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        samtools_container
    shell:
        'samtools view '
        '-T {input.fasta} '
        '-u '
        '{input.sam} '
        '2> {log} '
        '| '
        'samtools sort '
        '-l 0 '
        '-m 40G '
        '-O BAM '
        '--threads {threads} '
        '- '
        '> {output.bam} '
        '2>> {log} '

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
        '-a '
        '-x map-ont '
        '-t {threads} '
        '{input.fasta} '
        '{input.fq} '
        '> {output} '
        '2> {log}'


# chunk the assembly
rule list_contigs:
    input:
        'output/015_genome-chunks/chunk_{chunk}.fasta'
    output:
        temp('output/015_genome-chunks/chunk_{chunk}_contigs.txt')
    threads:
        1
    singularity:
        bbmap_container
    shell:
        'grep "^>" {input} | cut -d " " -f1 | sed -e \'s/>//g\' > {output}'

rule partition:
    input:
        'output/010_flye/assembly.fasta'
    output:
        temp(expand('output/015_genome-chunks/chunk_{chunk}.fasta',
                    chunk=all_chunks))
    params:
        outfile = 'output/015_genome-chunks/chunk_%.fasta',
        ways = n_chunks
    log:
        'output/015_genome-chunks/partition.log'
    threads:
        1
    singularity:
        bbmap_container
    shell:
        'partition.sh '
        '-Xmx800g '
        'in={input} '
        'out={params.outfile} '
        'ways={params.ways} '
        '2> {log}'


# assemble
rule flye_polish:
    input:
        fq = 'output/005_trim/ont_trimmed.fq',
        fa = 'output/010_flye/30-contigger/contigs.fasta'
    output:
        'output/020_flye-polish/polished_1.fasta'
    params:
        outdir = 'output/020_flye-polish',
        size = '10g'
    threads:
        min(128, multiprocessing.cpu_count())
    log:
        'output/logs/020_flye-polish.log'
    singularity:
        flye_container
    shell:
        'flye '
        '--nano-raw {input.fq} '
        '--polish-target {input.fa} '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&>> {log}'


rule flye:
    input:
        fq = 'output/005_trim/ont_trimmed.fq'
    output:
        # 'output/010_flye/assembly.fasta'
        'output/010_flye/30-contigger/contigs.fasta'
    params:
        outdir = 'output/010_flye',
        size = '10g'
    threads:
        min(128, multiprocessing.cpu_count())
    log:
        'output/logs/010_flye.log'
    singularity:
        flye_container
    shell:
        'flye '
        '--resume '
        '--nano-raw {input.fq} '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--stop-after repeat '
        '--threads {threads} '
        '&>> {log}'


rule gather_trimmed_reads:
    input:
        expand('output/005_trim/trimmed-chunk_{chunk}.fq',
               chunk=[str(x) for x in range(0, n_chunks)])
    output:
        'output/005_trim/ont_trimmed.fq'
    singularity:
        flye_container
    shell:
        'cat {input} > {output}'


rule remove_ont_adaptors:
    input:
        'output/005_trim/raw-chunk_{chunk}.fq'
    output:
        'output/005_trim/trimmed-chunk_{chunk}.fq'
    log:
        'output/logs/remove_ont_adaptors.{chunk}.log'
    threads:
        max(1, multiprocessing.cpu_count() // n_chunks)
    singularity:
        porechop
    shell:
        'porechop '
        '-i {input} '
        '-o {output} '
        '--verbosity 1 '
        '--threads {threads} '
        '--check_reads 1000 '
        '--discard_middle '
        '&> {log}'




rule chunk_raw_reads:
    input:
        ont_raw
    output:
        temp(expand('output/005_trim/raw-chunk_{chunk}.fq',
                    chunk=[str(x) for x in range(0, n_chunks)]))
    params:
        outfile = 'output/005_trim/raw-chunk_%.fq',
        ways = n_chunks
    log:
        'output/chunk_raw_reads.log'
    threads:
        1
    singularity:
        bbmap_container
    shell:
        'partition.sh '
        '-Xmx800g '
        'qin=33 '
        'in={input} '
        'out={params.outfile} '
        'ways={params.ways} '
        '2> {log}'



# generic BUSCO rule
rule busco:
    input:
        unpack(busco_resolver),
        lineage = 'data/metazoa_odb10'
    output:
        ('output/099_busco/'
         '{name}/run_metazoa_odb10/'
         'full_table.tsv'),
    log:
        Path(('output/logs/'
              'busco.{name}.log')).resolve()
    params:
        wd = 'output/099_busco',
        name = '{name}',
        fasta = lambda wildcards, input: Path(input.fasta).resolve(),
        lineage = lambda wildcards, input:
            Path(input.lineage).resolve()
    threads:
        multiprocessing.cpu_count()
    singularity:
        busco
    shell:
        'cd {params.wd} || exit 1 ; '
        'busco '
        '--force '
        '--in {params.fasta} '
        '--out {params.name} '
        '--lineage_dataset {params.lineage} '
        '--cpu {threads} '
        '--mode genome '
        '&> {log}'



# rule busco:
#     input:
#         unpack(busco_resolver),
#         lineage = 'data/busco/{lineage}'
#     output:
#         'output/099_busco/{lineage}/run_{name}/full_table_{name}.tsv'
#     log:
#         resolve_path('output/logs/099_busco/{lineage}-{name}.log')
#     params:
#         wd = 'output/099_busco/{lineage}',
#         fasta = lambda wildcards, input: resolve_path(input.fasta),
#         lineage = lambda wildcards, input: resolve_path(input.lineage),
#         tmpdir = tempfile.mkdtemp()
#     threads:
#         multiprocessing.cpu_count()
#     priority:
#         1
#     singularity:
#         busco_container
#     shell:
#         'cd {params.wd} || exit 1 ; '
#         'run_BUSCO.py '
#         '--force '
#         '--tmp_path {params.tmpdir} '
#         '--in {params.fasta} '
#         '--out {wildcards.name} '
#         '--lineage {params.lineage} '
#         '--cpu {threads} '
#         '--mode genome '
#         '&> {log}'

