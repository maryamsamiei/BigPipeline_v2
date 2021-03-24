# BigPipeline

This is the integration of EA-ML, EA-Wavelet, EPIMUTESTR and EA-Pathways (Reactomes and STRING Communities) pipelines. The output of each pipeline is a ranked list of genes based on the FDR values. 

## Installation
1. git clone https://github.com/LichtargeLab/BigPipeline.git
2. conda env create -f ./BigPipeline/environment.yml
3. conda install -n pyBigPipeline openjdk
4. conda activate pyBigPipeline


## Usage
Required arguments:
| Argument                | Descripion |
| ---------------------- |--------------------- |
| --VCF                | Path to annotated VCF file |
| --samples            |Path to two-column CSV file with sample IDs in the first column and patient labels (cases=1, controls=0) in the second column. There should be no header row in the csv|
| --savepath           | Path for output files |
| --cores              | number of cpus to use |

Optional arguments:
| Argument                 | Descripion |
| ---------------------- |--------------------- |
| --maxaf  | sets maximum allele frequency threshold (default: 1) |
| --minaf      |sets minimum allele frequency threshold (default: 0)|
| --transcript           |how to parse EA scores from different transcripts (options: canonical, max, mean, all / default: canonical)|
| --ref      | genome reference (options: hg19, hg38 / default: hg38) |
| --AC      | allele count (AC) threshold indicating which alleles from the reference population should be included in EA-Pathways output variant files / default: 5 |
| --refPop           | Path to text file containing reference population variants used to filter cohort VCF for EA-Pathways analysis (see "Additional Notes" below to see how to build it. (default: UKB_200K_WES.AC.AF.12102020.txt which is related to 200k UKB) |
| --pipeline           | which pipeline you want to run. (options: All, ML, Pathways, EAML, EPI, Wavelet/ default: All)|


## Command line example
```bash
python BigPipeline.py --VCF Path/to/vcf_file.vcf.gz --samples Path/to/samples_file.csv --savepath save/directory/ --cores 20 --maxaf 0.01 --AC 5 --pipeline ML
