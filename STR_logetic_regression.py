#!/usr/bin/env python
"""
Perform STR association tests using logistic regression.
"""

# Constants
MIN_STR_LENGTH = 6  # Minimum reference length for an STR

# Imports
import sys
import os
import warnings
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import vcf

from typing import List, Tuple, Optional, Dict

# Suppress warnings
warnings.filterwarnings("ignore")

# Add utility path
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "utils"))

import strtools.utils.common as common

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perform STR association tests using logistic regression (default).")

    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Output prefix", type=str, required=True)
    inout_group.add_argument("--fam", help="FAM file with phenotype info", type=str, required=True)
    inout_group.add_argument("--samples", help="File with list of samples to include", type=str)
    inout_group.add_argument("--exclude-samples", help="File with list of samples to exclude", type=str)

    pheno_group = parser.add_argument_group("Phenotypes")
    pheno_group.add_argument("--pheno", help="Phenotypes file (to use instead of --fam)", type=str)
    pheno_group.add_argument("--mpheno", help="Use (n+2)th column from --pheno", type=int, default=1)
    pheno_group.add_argument("--missing-phenotype", help="Missing phenotype code", type=str, default="-9")

    covar_group = parser.add_argument_group("Covariates")
    covar_group.add_argument("--covar", help="Covariates file", type=str)
    covar_group.add_argument("--covar-name", help="Names of covariates to load. Comma-separated", type=str)
    covar_group.add_argument("--covar-number", help="Column number of covariates to load. Comma-separated", type=str)
    covar_group.add_argument("--sex", help="Include sex from fam file as covariate", action="store_true")
    covar_group.add_argument("--cohort-pgc", help="Use cohort from PGC FIDs as a covariate", action="store_true")

    assoc_group = parser.add_argument_group("Association testing")
    assoc_group.add_argument("--region", help="Only process this region (chrom:start-end)", type=str)
    assoc_group.add_argument("--minmaf", help="Ignore bi-allelic sites with low MAF.", type=float, default=0.01)
    assoc_group.add_argument("--str-only", help="Only analyze STRs", action="store_true", default=True)
    assoc_group.add_argument("--remove-rare-str-alleles", help="Remove genotypes with alleles less than this freq", default=0.0, type=float)
    assoc_group.add_argument("--max-iter", help="Maximum number of iterations for logistic regression", default=100, type=int)

    fm_group = parser.add_argument_group("Fine mapping")
    fm_group.add_argument("--condition", help="Condition on this position chrom:start", type=str)

    return parser.parse_args()

def load_condition(vcf_file: str, condition: str, sample_order: List[str]) -> Tuple[List[int], List[str]]:
    """Load genotype data for a specific STR condition."""
    reader = vcf.Reader(open(vcf_file, "rb"))
    chrom, start = condition.split(":")
    region = f"{chrom}:{start}-{int(start) + 1}"
    for record in reader.fetch(region):
        if record.start == int(start):
            return load_genotypes(record, sample_order)
    common.ERROR("Could not find STR to condition on")

def get_assoc_type(alt_len: int = -1, name: Optional[str] = None) -> str:
    """Return string describing STR association type."""
    if alt_len >= 0:
        return f"{name}-length-{alt_len}"
    return name or "STR"

def print_header(outf, comment_lines: Optional[List[str]] = None):
    """Print header info for association output."""
    header = ["CHROM", "BP", "STR", "P", "OR", "SE", "CI95", "MAF", "NMISS", "ALLELE1", "ALLELE2"]
    outf.write("\t".join(header) + "\n")
    if comment_lines:
        for line in comment_lines:
            outf.write(f"# {line}\n")
    outf.flush()

def output_assoc(chrom: str, start: int, assoc: Dict[str, str], outf, assoc_type: str = "STR"):
    """Write association output."""
    if assoc is None:
        return
    items = [chrom, start, assoc_type, assoc["pval"], assoc["coef"], assoc["stderr"], assoc["CI"], assoc["maf"], assoc["N"], "CAG", "CAGCAG"]
    outf.write("\t".join(map(str, items)) + "\n")
    outf.flush()

def perform_association(data: pd.DataFrame, covarcols: List[str], maf: float = 1.0, 
                        exclude_samples: List[str] = [], maxiter: int = 100) -> Optional[Dict[str, str]]:
    """Perform STR association tests using logistic regression."""
    data = data[~data["sample"].isin(exclude_samples)]

    if data.empty:
        return None

    assoc = {"maf": f"{maf:.3f}", "N": str(data.shape[0])}
    
    try:
        model = sm.Logit(data["phenotype"], data[["intercept", "GT"] + covarcols])
        result = model.fit(disp=0, maxiter=maxiter)
        assoc["coef"] = f"{np.exp(result.params['GT']):.3f}"
        assoc["pval"] = f"{result.pvalues['GT']:.2E}"
        assoc["stderr"] = f"{np.exp(result.params['GT']) * result.bse['GT']:.3f}"
        assoc["CI"] = "-".join([f"{np.exp(ci):.3f}" for ci in result.conf_int().loc["GT"]])
    except Exception as e:
        common.MSG(f"Error in regression: {str(e)}")
        assoc.update({"coef": "NA", "pval": "NA", "stderr": "NA", "CI": "NA"})

    return assoc

def load_genotypes(record: vcf.model._Record, sample_order: List[str], 
                   use_alt_length: int = -1, rmrare: float = 0.0) -> Tuple[List[int], List[str]]:
    """Load STR genotypes from a record and return values in the sample order."""
    gtdata = {}
    exclude_samples = []
    alleles = [record.REF] + record.ALT
    afreqs = [(1 - sum(record.aaf))] + record.aaf

    for sample in record.samples:
        try:
            if use_alt_length > -1:
                gtdata[sample.sample] = sum(int(len(alleles[int(item)]) == use_alt_length) for item in sample.gt_alleles)
            else:
                if any(af < rmrare for af in [afreqs[int(item)] for item in sample.gt_alleles]):
                    exclude_samples.append(sample.sample)
                    gtdata[sample.sample] = sum(len(record.REF) for item in sample.gt_alleles)
                else:
                    gtdata[sample.sample] = sum(len(alleles[int(item)]) for item in sample.gt_alleles)
            if sample.sample not in gtdata:
                exclude_samples.append(sample.sample)
                gtdata[sample.sample] = 0
        except Exception as e:
            common.MSG(f"Error processing genotype: {str(e)}")
            exclude_samples.append(sample.sample)
            gtdata[sample.sample] = 0

    genotypes = [gtdata[s] for s in sample_order]
    return genotypes, exclude_samples

def load_pheno_data(fname: str, fam: bool = True, missing: str = "-9", mpheno: int = 1, sex: bool = False) -> pd.DataFrame:
    """Load phenotype data from FAM or pheno file."""
    if fam:
        columns = ["FID", "IID", "Father_ID", "Mother_ID", "sex", "phenotype"]
        data = pd.read_csv(fname, delim_whitespace=True, names=columns, dtype={"FID": str, "IID": str})
        if sex:
            data = data[data["sex"] != 0]  # 1=male, 2=female, 0=unknown
    else:
        columns = ["FID", "IID", "phenotype"]
        data = pd.read_csv(fname, delim_whitespace=True, names=columns, usecols=[0, 1, mpheno + 1])
        data = data[data["phenotype"] != missing]
        data["phenotype"] = data["phenotype"].apply(int) - 1  # Convert to 0/1 for logistic regression
    return data

def main():
    args = parse_arguments()

    # Load phenotype information
    common.MSG("Loading phenotype information...")
    if args.fam is not None:
        pdata = load_pheno_data(args.fam, fam=True, missing=args.missing_phenotype, sex=args.sex)
    elif args.pheno is not None:
        pdata = load_pheno_data(args.pheno, fam=False, missing=args.missing_phenotype, mpheno=args.mpheno)
    else:
        common.ERROR("Must specify phenotype using either --fam or --pheno")
    common.MSG(f"Loaded {pdata.shape[0]} samples...")

    # Load covariate information
    common.MSG("Loading covariate information...")
    covarcols = []
    if args.covar is not None:
        pdata, covarcols = add_covariates(pdata, args.covar, args.covar_name, args.covar_number)
    if args.sex:
        covarcols.append("sex")
    if args.cohort_pgc:
        pdata["cohort"] = pdata["FID"].apply(lambda x: x.split("*")[0])
        covarcols.append("cohort")
    common.MSG(f"Loaded {pdata.shape[0]} samples...")

    # Include/exclude samples
    common.MSG("Loading sample information...")
    if args.samples is not None:
        pdata = restrict_samples(pdata, args.samples, include=True)
    if args.exclude_samples is not None:
        pdata = restrict_samples(pdata, args.exclude_samples, include=False)
    common.MSG(f"Left with {pdata.shape[0]} samples...")

    # Setup VCF reader
    common.MSG("Set up VCF reader")
    reader = vcf.Reader(open(args.vcf, "rb"))

    # Set sample ID to FID_IID to match vcf
    common.MSG("Set up sample info")
    pdata["sample"] = pdata.apply(lambda x: f"{x['FID']}_{x['IID']}", axis=1)
    reader_samples = set(reader.samples)
    pdata = pdata[pdata["sample"].isin(reader_samples)]
    sample_order = list(pdata["sample"])
    pdata = pdata[["phenotype", "sample"] + covarcols]
    common.MSG(f"Left with {pdata.shape[0]} samples...")

    # Get data to condition on
    if args.condition is not None:
        cond_gt = load_condition(args.vcf, args.condition, sample_order)
        pdata["condition"] = cond_gt[0]
        covarcols.append("condition")

    # Prepare output file
    outf = sys.stdout if args.out == "stdout" else open(args.out, "w")
    print_header(outf, comment_lines=[" ".join(sys.argv)])

    # Perform association test for each record
    common.MSG(f"Performing associations... with covariates {str(covarcols)}")
    if args.region:
        reader = reader.fetch(args.region)
    for record in reader:
        # Check MAF
        aaf = sum(record.aaf)
        aaf = min([aaf, 1 - aaf])
        if aaf < args.minmaf:
            continue
        # Skip non-STRs
        if len(record.REF) < MIN_STR_LENGTH:
            continue
        # Extract genotypes in sample order, perform regression, and output
        common.MSG("   Loading genotypes...")
        gts, exclude_samples = load_genotypes(record, sample_order, rmrare=args.remove_rare_str_alleles)
        pdata["GT"] = gts
        pdata["intercept"] = 1
        common.MSG("   Performing association...")
        assoc = perform_association(pdata, covarcols, maf=aaf, exclude_samples=exclude_samples, maxiter=args.max_iter)
        common.MSG("   Outputting association...")
        output_assoc(record.CHROM, record.POS, assoc, outf, assoc_type=get_assoc_type(name=record.ID))

if __name__ == "__main__":
    main()

