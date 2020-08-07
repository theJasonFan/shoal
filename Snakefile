from os.path import join

### DATA
DATA_DIR = 'data'

# Experiments
url_path = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161{i}/DRR0161{i}_{LR}.fastq.gz'
fastq_fmt = join(DATA_DIR, 'DRR0161{i}', 'DRR0161{i}_{LR}.fastq.gz')
idxs = range(25,41)
#idxs = [25]
LR = [1,2]

# TX
tx_url = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz'
tx_path = join(DATA_DIR, 'athal.fa.gz')

### Quantify
EXP_NAME = config['exp_name']
OUTPUT_DIR = join('output', EXP_NAME)
READS_DIR = config['reads_dir']
SAMPLES_NAMES = config['sample_names']

TX_PATH = config['tx_path']
TX_INDEX = join(OUTPUT_DIR, '{}_index'.format(EXP_NAME))

QUANTS_DIR = join(OUTPUT_DIR, 'quants')
reads_1_fmt = join(READS_DIR, '{sample_name}', '{sample_name}_1.fastq.gz')
reads_2_fmt = join(READS_DIR, '{sample_name}', '{sample_name}_2.fastq.gz')

quant_dir_fmt = join(QUANTS_DIR, '{sample_name}')
eq_classes_gz_fmt = join(quant_dir_fmt, 'aux_info', 'eq_classes.txt.gz')
eq_classes_fmt = eq_classes_gz_fmt[:-3]

SHOAL = join('src-cpp', 'shoal')
SHOAL_OUTPUT_DIR = join(OUTPUT_DIR, 'shoal')
SHOAL_PRIOR = join(SHOAL_OUTPUT_DIR, 'prior.tsv')
shoal_quant_sf_fmt = join(SHOAL_OUTPUT_DIR, '{sample_name}_adapt.sf')

### SHOAL

rule all:
    input:
        expand(shoal_quant_sf_fmt, sample_name=SAMPLES_NAMES)

rule shoal_quant:
    input:
        prior = SHOAL_PRIOR,
         # TODO: figure out why snakmake won't let me write the following:
        # salmon_quant = quant_dir_fmt, 
    output:
        shoal_quant_sf_fmt
    params:
        sf=quant_dir_fmt # hack for comment above
    shell:
        '''
        ./{SHOAL} -p {input.prior} -s {params.sf} -o {output} -t adapt-prior
        '''

SHOAL_PY = join('src-py', 'shoal.py')
rule shoal_prior:
    input:
        expand(eq_classes_fmt, sample_name=SAMPLES_NAMES)
    output:
        prior=SHOAL_PRIOR,
    params:
        samples=','.join(SAMPLES_NAMES),
        basepath=QUANTS_DIR,
        outdir=directory(SHOAL_OUTPUT_DIR)
    shell:
        '''
        python {SHOAL_PY} create-prior\
            --samples {params.samples} \
            --basepath {params.basepath} \
            --outdir {params.outdir}
        '''


rule unzip_all:
    input:
        expand(eq_classes_fmt, sample_name=SAMPLES_NAMES)

rule unzip:
    input:
        eq_classes_gz_fmt
    output:
        eq_classes_fmt
    shell:
        # Don't use --keep, the following is bacwards compatible with gzip < v1.6
        '''
        gunzip < {input} > {output}
        '''

rule quant_all:
    input:
        expand(eq_classes_gz_fmt, sample_name=SAMPLES_NAMES)
        
###
# SALMON params
# `--dumpEqWeights` to dump equivalence class file
# `--rangeFactorizationBins 1` to use traditional rich eq classes
rule quant:
    input:
        seq1=reads_1_fmt,
        seq2=reads_2_fmt,
        index=TX_INDEX,
    output:
        #outdir = directory(quant_dir_fmt),# possible bug w snakemake
        eq_classes_gz=eq_classes_gz_fmt,
    params:
        cores = config.get('salmon_cores', 8),
        outdir = directory(quant_dir_fmt), # sidestepping possible snakemake bug
    shell:
        '''
        salmon quant -i {input.index} -l A \
            -1 data/DRR016125/DRR016125_1.fastq.gz \
            -2 data/DRR016125/DRR016125_2.fastq.gz \
            -o {params.outdir} \
            --rangeFactorizationBins 1\
            -p {params.cores} \
            --validateMappings \
            --dumpEqWeights
        '''

rule index:
    input:
        TX_PATH
    output:
        directory(TX_INDEX)
    shell:
        '''
        salmon index -t {input} -i {output}
        '''

### Compile shoal
rule compile:
    input:
        makefle=join('src-cpp', 'Makefile'),
        d=('src-cpp')
    output:
        SHOAL
    shell:
        '''
        make --directory={input.d}
        '''

