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
| --refPop           | Path to text file containing reference population variants used to filter cohort VCF for EA-Pathways analysis (see "Additional Notes" below to see how to build it.)|
| --savepath           | Path for output files |
| --cores              | number of cpus to use |

Optional arguments:
| Argument                 | Descripion |
| ---------------------- |--------------------- |
| --maxaf  | sets maximum allele frequency threshold for ML pipelines (default: 0.01) |
| --minaf      |sets minimum allele frequency threshold for ML pipelines (default: 0)|
| --minAC      | Minimum allele count (AC) cutoff based on the reference population for EA-Pathways analysis / default: 1 |
| --maxAC      | Maximum allele count (AC) threshold based on the reference population for EA-Pathways analysis / default: 5 |
| --Ann      | Variant annotation pipeline used (options: ANNOVAR, VEP / default: VEP) |
| --transcript           |How to parse EA scores from different transcripts (options: canonical, max, mean, all / default: canonical)|
| --ref      | Genome reference (options: hg19, hg38 / default: hg38) |
| --pipeline           | Which pipeline you want to run. (options: All, ML, Reactome, STRING, GoTerms, sigma, EAML, EPI, Wavelet/ default: All)|
| --chrX       | Whether to include (1) or exclude (0) sex chromosomes in analysis (options: 1, 0 / default: 1 )|
| --JVMmemory  | Memory needed for each Weka JVM. (use Xmx3g if you got "BlockingIOError: Resource temporarily unavailable" / default: Xmx2g)|
| --GeneLength  | gene length file path/default is './refs/gene_length.csv'|

Note: minaf and maxaf set thresholds for ML pipelines (EA-ML, EPIMUTESTR and EA-Wavelet) while minAC and maxAC set thresholds for EA-Pathways pipelines (Reactomes and STRING Communities). 

## Command line example
```bash
#set your working directory to BigPipeline
cd ./BigPipeline
#run Big_Pipeline.py
python ./src/Big_Pipeline.py --VCF Path/to/vcf_file.vcf.gz --samples Path/to/samples_file.csv --savepath save/directory/ --cores 20 --maxaf 0.01 --minAC 3 --maxAC 7 --pipeline All
```

## Additional Notes
* Reference population variant file format:
If reference population variants are in VCF format, then user can use the following command to generate the reference population variant file:
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AF\n' reference_population.vcf.gz > reference_populaton_INFO.txt.
Otherwise, reference txt file should have the following columns for each variant and no header row:
  * Column 1, "CHROM": variant chromosome location 
  * Column 2, "POS": variant start location 
  * Column 3, "REF": variant reference allele 
  * Column 4, "ALT": variant alternative allele 
  * Column 5, "AC": variant allele count in reference population 
  * Column 6, "AF": variant allele frequency in reference population 

* Weka must be downloaded separately. Click [here](https://waikato.github.io/weka-wiki/downloading_weka/) to install Weka 3.8.0+. 
Put Weka folder in the BigPipeline directory.
* For a sample size around 7000, 400GB memory is needed. It's highly recommended to run BigPipeline using at least 20 cpus in parallel.

