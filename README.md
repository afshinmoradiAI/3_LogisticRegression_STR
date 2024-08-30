STR Association Test Script
This repository contains a Python script designed to perform Short Tandem Repeat (STR) association tests using logistic regression. The script is optimized for handling genotype data in VCF format and allows for the incorporation of phenotype and covariate data. It also supports conditioning on specific genomic positions and performing allele-based and allele-length-based tests.
Features
• VCF File Loading: The script loads VCF files and processes STR data.
• Phenotype Data Handling: The script can load and process phenotype data from FAM or custom phenotype files.
• Covariate Integration: Covariates can be added to the phenotype data to control for confounding factors.
• Sample Filtering: Allows for inclusion or exclusion of specific samples based on an input file.
• Conditioning: Supports conditioning on specific genomic positions to account for known associations.
• Logistic Regression Association Testing: Performs logistic regression for association testing, with support for various allele-based and allele-length-based tests.
• Output: Results are outputted in a tab-delimited format, suitable for downstream analysis.
Dependencies
The script requires the following Python libraries:

- numpy
- pandas
- statsmodels
- vcf (PyVCF)
- argparse
- logging
Installation
To install the necessary dependencies, you can use the following command:
pip install numpy pandas statsmodels PyVCF argparse
Usage
To run the script, use the following command in your terminal:
python str_association.py --vcf <path_to_vcf_file> --out <output_file_prefix> --fam <fam_file> --covar <covar_file>
Command-Line Arguments
--vcf: Path to the input VCF file.
--out: Prefix for the output file.
--fam: FAM file containing phenotype information.
--samples: File with a list of samples to include.
--pheno: Custom phenotype file (to use instead of FAM).
--covar: File containing covariate data.
--covar-name: Names of covariates to include (comma-separated).
--covar-number: Column numbers of covariates to include (comma-separated).
--sex: Include sex as a covariate.
--region: Process only this specific region (chrom:start-end).
--minmaf: Minimum Minor Allele Frequency (MAF) for biallelic sites.
--condition: Condition on this position (chrom:start).
Example
python str_association.py --vcf example.vcf --out results --fam example.fam --covar example.covar --region chr1:100000-200000
Output
The script outputs a tab-delimited file containing the results of the association tests. The columns include:

- CHROM: Chromosome
- BP: Base pair position
- SNP: SNP identifier
- P: p-value from the logistic regression
- OR: Odds ratio
- SE: Standard error
- CI95: 95% confidence interval
- MAF: Minor allele frequency
- NMISS: Number of missing samples
- ALLELE1: Reference allele
- ALLELE2: Alternative allele
Logging
The script uses Python's built-in logging module to provide information about the processing steps. The logs include data loading, filtering, and association testing steps.
License
This project is licensed under the MIT License.
Contact
For questions or issues, please contact [your email] or create an issue on this repository.
