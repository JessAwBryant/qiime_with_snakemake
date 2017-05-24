# qiime_with_snakemake

Author: Jessica A Bryant (with substantial help from John Eppley)
Affiliation: UH and MIT

Aim: A simple snakemake workflow to process MiSeq paired-end SSU amplicon data using QIIME.

Required Files:

	config.yaml: 
		This file tells snakemake where the data files are. Two files and corresponding file paths are required:
		
		1) A single directory containing all raw paired-end MiSeq files (forward & reverse reads plus index files).
		2) A Qiime mapping file.
                   Formatting information Available at: http://qiime.org/documentation/file_formats.html

    alternative.yaml:
    	 This is an alternative configuration file to set up an environment to run Qiime rules. 
    	 This is necessary because Snakemake runs in python3 but Qiime currently only runs in Python2.7.             
                 
Date: March 22, 2017
    To run in my directory:
    
     1. Activate a python3 environment. (for more info: https://conda.io/docs/py2or3.html)
        ex: source activate py35
        
     2. Edit the config.yml with the correct mapping file path and
         the directory containing the raw MiSeq data.
         
     3. place this Snakefile, updated config.yml and alternative.yml in the same directory. 
     Then from this directory type: snakemake --use-conda
