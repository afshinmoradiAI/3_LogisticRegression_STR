# STR Association Test Script

## Introduction
- This script performs logistic regression-based association tests on STRs within a VCF file, considering various covariates and phenotype data. The results include odds ratios, p-values, confidence intervals, and more, making it a valuable tool for genetic association studies.
  
## Features
- **VCF File Parsing:** Reads genetic variant data from VCF files.
- **Phenotype Data Loading:** Supports loading phenotype data from FAM or custom phenotype files.
- **Covariate Inclusion:** Allows adding covariates, including sex and cohort data, to the analysis.
- **Sample Filtering:** Includes options to restrict or exclude specific samples from the analysis.
- **Conditional Analysis:** Supports conditional logistic regression based on specific genetic variants.
- **Multiple Test Types:** Performs allele-based and allele-length-based association tests.

## Installation
To run the script, you need to have Python 3.x installed, along with the required dependencies.

1. **Clone the repository:**
   ```bash
   git clone https://github.com/your-username/str-association-test.git
   cd str-association-test
   ```

2. **Install the required Python packages:**
   ```bash
   pip install -r requirements.txt
   ```
   Install the dependencies via `pip`:

   ```bash
   pip install numpy pandas statsmodels pyVCF argparse logging
   ```
 
Usage
To use the script, run it from the command line with the appropriate arguments:

  ```bash
  python str_association_test.py --vcf input.vcf --fam input.fam --out results.txt --covar covariates.txt --region 1:10000-20000
  ```

## Example Command

```bash
python str_association_test.py --vcf data/variants.vcf --fam data/sample.fam --out results/association_results.txt --covar data/covariates.txt --region 1:150000-250000 --condition 1:160000 --minmaf 0.05
```
Input Files
- **VCF File:** Contains genetic variants.
- **FAM File:** Includes phenotype information.
- **Phenotype File:** (Optional) An alternative to FAM files for phenotype data.
- **Covariates File:** Contains additional covariate data, such as sex or cohort information.
Output
The script outputs a tab-delimited text file containing the results of the association tests. The key columns include:

- `CHROM`: Chromosome number
- `BP`: Base pair position
- `SNP`: Identifier for the STR
- `P`: p-value for the association test
- `OR`: Odds ratio
- `SE`: Standard error
- `CI95`: 95% confidence interval
- `MAF`: Minor allele frequency
- `NMISS`: Number of samples included
- `ALLELE1`, `ALLELE2`: Reference and alternate alleles
Options
- `--vcf`: Path to the input VCF file (required).
- `--out`: Output file prefix or path (required).
- `--fam`: Path to the FAM file with phenotype information.
- `--pheno`: Path to the phenotype file (alternative to FAM file).
- `--covar`: Path to the covariates file.
- `--region`: Specific genomic region to analyze (e.g., `chr:start-end`).
- `--minmaf`: Minimum minor allele frequency to consider.
- `--condition`: Chromosome position to condition on for logistic regression.
- `--samples`: File with a list of samples to include.


License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
