from os.path import join

### DATA
DATA_DIR = 'data'

# Experiments
url_path = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161{i}/DRR0161{i}_{LR}.fastq.gz'
fastq_fmt = join(DATA_DIR, 'DRR0161{i}', 'DRR0161{i}_{LR}.fastq.gz')
idxs = range(25,41)
idxs = [25]
LR = [1,2]

# TX
tx_url = 'ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz'
tx_path = join(DATA_DIR, 'athal.fa.gz')

### Quantify
OUTPUT_DIR = 'output'
QUANTS_DIR = join(OUTPUT_DIR, 'quants')
quant_dir_fmt = join(QUANTS_DIR, 'DRR0161{i}')
eq_classes_gz_fmt = join(quant_dir_fmt, 'aux_info', 'eq_classes.txt.gz')
eq_classes_fmt = eq_classes_gz_fmt[:-3]


athal_index = join(OUTPUT_DIR, 'athal_index')



### SHOAL

rule unzip_all:
    input:
        expand(eq_classes_fmt, i=idxs)

rule unzip:
    input:
        eq_classes_gz_fmt
    output:
        eq_classes_fmt
    shell:
        # Don't use --keep, this is bacwards compatible with gzip < v1.6
        '''
        gunzip < {input} > {output}
        '''

rule quant_all:
    input:
        expand(eq_classes_gz_fmt, i=idxs)
        
###
# SALMON params
# `--dumpEqWeights` to dump equivalence class file
# `--rangeFactorizationBins 1` to use traditional rich eq classes
rule quant:
    input:
        # Dirty trick for incremental string fmt
        seq1=fastq_fmt.replace('{LR}', '1'),
        seq2=fastq_fmt.replace('{LR}', '2'),
        index=athal_index,
    output:
        outdir = directory(quant_dir_fmt),
        eq_classes_gz=eq_classes_gz_fmt,
    params:
        cores = config.get('salmon_cores', 8)
    shell:
        '''
        salmon quant -i {input.index} -l A \
            -1 data/DRR016125/DRR016125_1.fastq.gz \
            -2 data/DRR016125/DRR016125_2.fastq.gz \
            -o {output.outdir} \
            --rangeFactorizationBins 1\
            -p {params.cores} \
            --validateMappings \
            --dumpEqWeights
        '''

rule index:
    input:
        tx_path
    output:
        directory(athal_index)
    shell:
        '''
        salmon index -t {input} -i {output}
        '''
        
### Download data from Salmon tutorial
rule download_all:
    input:
        expand(fastq_fmt, i=idxs, LR=LR),
        tx_path

rule download:
    output:
        fastq_fmt
    params:
        url = lambda w: url_path.format(i=w.i, LR=w.LR)
    shell:
        '''
        wget {params.url} -O {output}
        '''

rule download_tx:
    output:
        tx_path
    params:
        url = tx_url
    shell:
        '''
        wget  {params.url} -O {output}
        '''
