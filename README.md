How to Run the STR Logistic Regression Script

To run the STR logistic regression script, follow the steps below:

### Prerequisites
Ensure you have Python installed and the necessary Python packages. You can install the required packages using pip:

```bash
pip install numpy pandas statsmodels PyVCF
```

### Running the Script
Use the following command to run the script:

```bash
python your_script_name.py --vcf example.vcf --out result --fam example.fam --minmaf 0.05 --str-only --max-iter 200

```
In this example:

    --vcf example.vcf specifies the input VCF file.
    --out result sets the prefix for the output files.
    --fam example.fam specifies the FAM file with phenotype information.
    --minmaf 0.05 sets the minimum minor allele frequency (MAF) for filtering (optional).
    --str-only ensures that only STRs are analyzed.
    --max-iter 200 sets the maximum number of iterations for logistic regression (optional).

    
Output Explanation

### Output Files
The script generates several output files with a prefix specified by the `--out` argument. These files might include:

1. **Summary Statistics File** (`<output_prefix>_summary.txt` or similar):
   - This file contains the summary statistics of the logistic regression tests performed for each STR. 
   - Typical columns might include:
     - **STR ID**: Identifier for the STR being tested.
     - **P-value**: Significance level of the association test.
     - **Odds Ratio (OR)**: The strength of association between the STR and the trait.
     - **Confidence Intervals (CI)**: The range within which the true odds ratio is expected to fall.
     - **Effect Size**: The estimated effect size of the STR on the trait.

2. **Detailed Results File** (`<output_prefix>_detailed.txt` or similar):
   - This file contains detailed logistic regression results for each STR, possibly including:
     - **Log-Likelihood**: A measure of how well the model fits the data.
     - **Coefficients**: The estimated coefficients for each predictor in the model.
     - **Standard Errors**: The standard error of each coefficient estimate.
     - **Z-scores**: The test statistic for each coefficient.
     - **P-values**: The significance level of each coefficient.

3. **Log File** (`<output_prefix>.log`):
   - A log file that tracks the progress of the script's execution. It might include information about the number of STRs tested, any warnings or errors encountered during the run, and summary messages about the execution time and resource usage.

### Understanding the Output

1. **P-value**:
   - A small p-value (typically ≤ 0.05) indicates strong evidence against the null hypothesis, suggesting that the STR is significantly associated with the trait.

2. **Odds Ratio (OR)**:
   - An OR greater than 1 suggests that the presence of a specific STR allele increases the likelihood of the trait.
   - An OR less than 1 suggests that the STR allele decreases the likelihood of the trait.

3. **Confidence Intervals (CI)**:
   - If the confidence interval for the OR does not include 1, the association is considered statistically significant.

### Example of Output Interpretation

If you ran the script with the following command:

```bash
python STR_logetic_regression.py --vcf data/sample.vcf --out results/output
```

You might get output files like:

- `results/output_summary.txt`
- `results/output_detailed.txt`
- `results/output.log`

In the `output_summary.txt` file, a row might look like this:

| STR_ID | P-value | OR   | CI_Lower | CI_Upper | Effect_Size |
|--------|---------|------|----------|----------|-------------|
| STR123 | 0.002   | 1.45 | 1.10     | 1.90     | 0.37        |

**Interpretation**:
- The STR `STR123` has a p-value of 0.002, indicating a statistically significant association with the trait.
- The odds ratio (OR) is 1.45, meaning the presence of this STR allele is associated with a 45% increase in the likelihood of the trait.
- The confidence interval (1.10 to 1.90) suggests that the true OR is likely to be within this range.
- The effect size is 0.37, quantifying the magnitude of the association.

