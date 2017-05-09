"""

Author: Jessica A Bryant (with substantial help from John E.)
Affiliation: UH and MIT
Aim: A simple snakemake workflow to process MiSeq paired-end SSU amplicon data.

Config file: 1) The workflow assumes all your samples came from one MiSeq run, so the input
             is the directory containing raw fastqs (forward & reverse reads plus index files).
             
              2) Qiime requires a mapping file
                 http://qiime.org/documentation/file_formats.html
                 
organizational aspirations: http://snakemake.readthedocs.io/en/latest/snakefiles/deployment.html#integrated-package-management                 

Date: March 22, 2017
    To run:
     1. Activate a python3 environment. (for more info: https://conda.io/docs/py2or3.html)
        ex: source activate py35
     2. Edit the config.yml with the correct mapping file path and
         the directory containing the raw MiSeq data.
     3. place this Snakefile, updated config.yml and alternative.yml in the same directory. 
     Then from this directory type: snakemake --use-conda

"""

configfile: "config.yaml"
SAMPLES, = glob_wildcards(config['raw_data_directory']+"/{reads}_R1_001.fastq")

for reads in SAMPLES:
  print("Sample[s] " + reads + " will be processed")

#declare the first target file as output
rule all: 
  input: 
      expand("{reads}.stats", reads=SAMPLES),
      expand("trimmed/{reads}_paired.fastq", reads=SAMPLES),
      expand("trimmed/{reads}_pear.assembled.fastq", reads=SAMPLES),
      expand("trimmed/{reads}_pear_reads_to_keep.txt", reads=SAMPLES),
      "trimmed/all_index.fastq",
      "trimmed/all_reads.fastq",
      "trimmed/seqs.uchime.fna",
      "trimmed/otus/otu_table_mc2_w_tax.biom",
      "trimmed/otus/status.txt"

rule get_sequencing_stats:
    input:
       fwd="folder/{reads}_R1_001.fastq",
       rev="folder/{reads}_R2_001.fastq"
    output:
        "{reads}.stats"
    conda:
         "/mnt/lysine/jbryant/dipole/snakemake/alternative.yaml"
    message: """Getting Stats"""
    shell:
        "count_seqs.py -i {input[0]} > {output}; count_seqs.py -i {input[1]} >> {output}"

rule trimming:
   input: 
       fwd="folder/{reads}_R1_001.fastq",
       rev="folder/{reads}_R2_001.fastq"
   output:
       fwd="trimmed/{reads}_paired.fastq",
       single="trimmed/{reads}_unpaired.fastq",
       rev="trimmed/{reads}_r_paired.fastq",
       single_rev="trimmed/{reads}_r_unpaired.fastq"
   message: """--- Trimming reads."""
   conda:
        "/mnt/lysine/jbryant/dipole/snakemake/alternative.yaml"
   shell:
       """
       trimmomatic PE -trimlog {wildcards.reads}_trimmed.log {input.fwd} {input.rev} {output.fwd} {output.single} {output.rev} {output.single_rev} ILLUMINACLIP:/home/jbryant/programs/Trimmomatic-0.36/TruSeq2-jb.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100
       """

rule assemble:
    input:
        fwd = "trimmed/{reads}_paired.fastq",
        rev = "trimmed/{reads}_r_paired.fastq"
    output: 
        asmba = "trimmed/{reads}_pear.assembled.fastq"
    params: 
        res="trimmed/{reads}_pear"
    conda:
        "/mnt/lysine/jbryant/dipole/snakemake/alternative.yaml"

    message: """--- Assembling reads"""
    log:
        "trimmed/{reads}.log"
    shell:
       """
       pear -f  {input.fwd}  -r {input.rev}  -o {params.res} > {log}
       """
       
rule barcodes_to_keep:
    input: 
        "trimmed/{reads}_pear.assembled.fastq"
    output:
        "trimmed/{reads}_pear_reads_to_keep.txt"
    message: """--- screen barcodes"""
    shell:
       """
       sed -n '1~4'p {input} | sed 's/^@//g' > {output}
       """

rule screen_index_fastq:
    input:
        fastq="folder/{reads}_I1_001.fastq",
        readlist="trimmed/{reads}_pear_reads_to_keep.txt"
    output:
        "trimmed/{reads}_screened_index.fastq"
    message: """--- cleaning up barcodes"""
    conda:
        "/mnt/lysine/jbryant/dipole/snakemake/alternative.yaml"
    shell:
       """ 
       filter_fasta.py -f {input.fastq} -o {output} -s {input.readlist}
       """

#if samples are spread across more than one file, this will concatinate them.
rule cat_samples:
    input:
        expand("trimmed/{reads}_pear.assembled.fastq",reads=SAMPLES)
    output:
        "trimmed/all_reads.fastq"
    message: 
        """--- concat sequences"""
    shell:
        """
        cat {input} > {output}
        """

#concatinate index files
rule cat_index:
    input:
        expand("trimmed/{reads}_screened_index.fastq",reads=SAMPLES)
    output:
        "trimmed/all_index.fastq"
    shell:
        """
        cat {input} > {output}
        """
        
rule split_samples:
    input:
        ins = "trimmed/all_reads.fastq",
        index = "trimmed/all_index.fastq",
        mapping = config['mapping_file']
    output:
        "trimmed/seqs.fna"
    params:
        folder = "trimmed/"
    message: 
        """--- splitting_samples"""
    conda:
        "/mnt/lysine/jbryant/dipole/snakemake/alternative.yaml"
    shell:
        """
        split_libraries_fastq.py -i {input.ins} --barcode_type 12 -m {input.mapping} -o {params.folder} -b  {input.index}
        touch {output}
        """
        
rule remove_chimeras:
    input:
        "trimmed/seqs.fna"
    output:
        "trimmed/seqs.uchime.fna"
    params:
        reference_database = "/home/jbryant/anaconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta"
    log:
        "trimmed/seqs.uchime.log"
    conda:
        "/mnt/lysine/jbryant/dipole/snakemake/alternative.yaml"
    shell:
       """
       vsearch --uchime_ref {input} --db {params.reference_database} --nonchimeras {output} -log {log}
       """
       
rule pick_open_reference_otus:
    input:
        "trimmed/seqs.uchime.fna"
    params:
        output_directory = "trimmed/otus",
        silva_params = "/mnt/lysine/jbryant/dipole/qiime_workflow/trimmed/silva-128-params.txt",
        silva_database = "/slipstream/home/jbryant/databases/SILVA_128_For_QIIME/SILVA_128_QIIME_release/rep_set/rep_set_all/97/97_otus.fasta"
    output:
        "trimmed/otus/otu_table_mc2_w_tax.biom"
    conda:
        "/mnt/lysine/jbryant/dipole/snakemake/alternative.yaml"
    shell:
        """
        pick_open_reference_otus.py -o {params.output_directory} -i {input} -p {params.silva_params} -r {params.silva_database} -n new  -f
        touch {output}
        """
        
rule basic_diversity_calculations:
    input:
        "trimmed/otus/otu_table_mc2_w_tax.biom"
    params:
        output_dir = "trimmed/otus/basic_div"
    output:
        "trimmed/otus/status.txt"
    conda:
        "/mnt/lysine/jbryant/dipole/snakemake/alternative.yaml"
    shell:
        """
        core_diversity_analyses.py -o {params.output_dir} -i {input} -m /mnt/lysine/jbryant/dipole/dipole_mapping.tsv -t /mnt/lysine/jbryant/dipole/snakemake/trimmed/otus/rep_set.tre -e 73741
        touch {output}
        """