#!/usr/bin/env python

MIN_STR_LENGTH = 6 

# Imports
import sys
import os
import warnings
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import vcf
import logging

# Custom imports - adjust based on your setup
import strtools.utils.common as common

# Suppress warnings
warnings.filterwarnings("ignore")

def setup_logging():
    """Setup logging for the script."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_vcf_file(vcf_file):
    """Load VCF file using the VCF reader."""
    return vcf.Reader(open(vcf_file, "rb"))

def load_pheno_data(fname, fam=True, missing=-9, mpheno=1, sex=False):
    """Load phenotype data from FAM or phenotype file."""
    if fam:
        data = pd.read_csv(fname, delim_whitespace=True, names=["FID", "IID", "Father_ID", "Mother_ID", "sex", "phenotype"])
        if sex:
            data = data[data["sex"] != 0]  # 1=male, 2=female, 0=unknown
    else:
        data = pd.read_csv(fname, delim_whitespace=True, names=["FID", "IID", "phenotype"], usecols=[0, 1, 1 + mpheno])
    
    data = data[data["phenotype"].apply(str) != missing]
    data["phenotype"] = data["phenotype"].apply(int) - 1  # Convert to 0/1
    return data

def add_covariates(data, covar_file, covar_name, covar_number):
    """Add covariates to phenotype data frame."""
    default_cols = ["FID", "IID"]
    if covar_name:
        colnames = default_cols + covar_name.split(",")
        cov = pd.read_csv(covar_file, delim_whitespace=True, usecols=colnames)
    elif covar_number:
        colnames = default_cols + ["C" + item for item in covar_number.split(",")]
        cov = pd.read_csv(covar_file, delim_whitespace=True, names=colnames, usecols=list(range(2)) + [1 + int(item) for item in covar_number.split(",")])
    else:
        cov = pd.read_csv(covar_file, delim_whitespace=True)
        if "FID" not in cov.columns:
            cov.columns = default_cols + cov.columns[len(default_cols):]
    
    data = pd.merge(data, cov, on=["FID", "IID"])
    covarcols = [col for col in data.columns if col not in default_cols + ["Father_ID", "Mother_ID", "sex", "phenotype"]]
    return data, covarcols

def restrict_samples(data, sample_file, include=True):
    """Include or exclude specific samples."""
    samples = pd.read_csv(sample_file, names=["FID", "IID"], delim_whitespace=True)
    samples = samples.applymap(str)
    if include:
        data = pd.merge(data, samples, on=["FID", "IID"])
    else:
        data = pd.merge(data, samples, on=["FID", "IID"], how="left", indicator=True)
        data = data[data["_merge"] == "left_only"]
        data.drop("_merge", axis=1, inplace=True)
    return data

def load_condition(vcf_file, condition, sample_order):
    """Load condition data from VCF."""
    reader = load_vcf_file(vcf_file)
    chrom, start = condition.split(":")
    region = f"{chrom}:{start}-{int(start) + 1}"
    reader.fetch(region)
    for record in reader:
        if record.start == int(start):
            return load_gt(record, sample_order, is_str=True)
    raise ValueError("Could not find SNP to condition on")

def load_gt(record, sample_order, is_str=True, use_alt_num=-1, use_alt_length=-1, rmrare=0):
    """Load genotypes from a record and return values in the sample order."""
    gtdata = {}
    exclude_samples = []
    alleles = [record.REF] + record.ALT
    afreqs = [(1 - sum(record.aaf))] + record.aaf
    
    for sample in record:
        try:
            if use_alt_num > -1:
                sum_alleles = sum([int(int(item) == use_alt_num) for item in sample.gt_alleles])
                gtdata[sample.sample] = sum_alleles
            elif use_alt_length > -1:
                sum_alleles_len = sum([int(len(alleles[int(item)]) == use_alt_length) for item in sample.gt_alleles])
                gtdata[sample.sample] = sum_alleles_len
            else:
                f1, f2 = [afreqs[int(item)] for item in sample.gt_alleles]
                if f1 < rmrare or f2 < rmrare:
                    exclude_samples.append(sample.sample)
                    gtdata[sample.sample] = sum([len(record.REF) for item in sample.gt_alleles])
                else:
                    gtdata[sample.sample] = sum([len(alleles[int(item)]) for item in sample.gt_alleles])
        except Exception as e:
            logging.warning(f"Error processing genotype for sample {sample.sample}: {e}")
            exclude_samples.append(sample.sample)
            gtdata[sample.sample] = 2 * len(record.REF)
    
    genotypes = [gtdata[s] for s in sample_order]
    return genotypes, exclude_samples

def perform_association(data, covarcols, maf=1.0, exclude_samples=[], maxiter=100):
    """Perform logistic regression association tests."""
    data = data[~data["sample"].isin(exclude_samples)]
    
    if data.empty:
        return None

    assoc = {}
    formula = "phenotype ~ GT+" + "+".join(covarcols)
    assoc["maf"] = f"{maf:.3f}"
    assoc["N"] = data.shape[0]
    
    try:
        logit_model = sm.Logit(data["phenotype"], data[["intercept", "GT"] + covarcols]).fit(disp=0, maxiter=maxiter)
        assoc["coef"] = f"{np.exp(logit_model.params['GT']):.3f}"
        assoc["pval"] = f"{logit_model.pvalues['GT']:.2E}"
        assoc["stderr"] = f"{np.exp(logit_model.params['GT']) * logit_model.bse['GT']:.3f}"
        assoc["CI"] = "-".join([f"{np.exp(logit_model.conf_int().loc['GT', ind]):.3f}" for ind in [0, 1]])
    except Exception as e:
        logging.error(f"Logistic regression failed: {e}")
        assoc.update({"coef": "NA", "pval": "NA", "stderr": "NA", "CI": "NA"})
    
    return assoc

def print_header(outf, comment_lines=[]):
    """Print header info for association output."""
    header = ["CHROM", "BP", "SNP", "P", "OR", "SE", "CI95", "MAF", "NMISS", "ALLELE1", "ALLELE2"]
    outf.write("\t".join(header) + "\n")
    outf.flush()

def output_assoc(chrom, start, assoc, outf, assoc_type="STR"):
    """Write association output."""
    if assoc is None:
        return
    items = [chrom, start, assoc_type, assoc["pval"], assoc["coef"], assoc["stderr"], assoc["CI"], assoc["maf"], assoc["N"], "CAG", "CAGCAG"]
    outf.write("\t".join(map(str, items)) + "\n")
    outf.flush()

def main():
    setup_logging()
    
    parser = argparse.ArgumentParser(description="Perform STR association tests using logistic regression")
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input VCF file", type=str)
    inout_group.add_argument("--out", help="Output prefix", type=str)
    inout_group.add_argument("--fam", help="FAM file with phenotype info", type=str)
    inout_group.add_argument("--samples", help="File with list of samples to include", type=str)
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
    assoc_group.add_argument("--allele-tests", help="Also perform allele-based tests using each separate allele", action="store_true")
    assoc_group.add_argument("--minmaf", help="Ignore bi-allelic sites with low MAF.", type=float, default=0.01)
    assoc_group.add_argument("--remove-rare-str-alleles", help="Remove genotypes with alleles less than this freq", default=0.0, type=float)
    assoc_group.add_argument("--max-iter", help="Maximum number of iterations for logistic regression", default=100, type=int)
    fm_group = parser.add_argument_group("Fine mapping")
    fm_group.add_argument("--condition", help="Condition on this position chrom:start", type=str)
    args = parser.parse_args()

    # Load phenotype information
    logging.info("Loading phenotype information...")
    pdata = load_pheno_data(args.fam if args.fam else args.pheno, fam=bool(args.fam), missing=args.missing_phenotype, sex=args.sex)
    logging.info(f"Loaded {pdata.shape[0]} samples")

    # Load covariate information
    logging.info("Loading covariate information...")
    covarcols = []
    if args.covar:
        pdata, covarcols = add_covariates(pdata, args.covar, args.covar_name, args.covar_number)
    if args.sex:
        covarcols.append("sex")
    if args.cohort_pgc:
        pdata["cohort"] = pdata["FID"].apply(lambda x: x.split("*")[0])
        covarcols.append("cohort")
    logging.info(f"Loaded {pdata.shape[0]} samples with covariates")

    # Include/exclude samples
    logging.info("Loading sample information...")
    if args.samples:
        pdata = restrict_samples(pdata, args.samples, include=True)
    logging.info(f"Left with {pdata.shape[0]} samples after filtering")

    # Setup VCF reader
    logging.info("Setting up VCF reader...")
    reader = load_vcf_file(args.vcf)

    # Set sample ID to FID_IID to match VCF
    pdata["sample"] = pdata.apply(lambda x: f"{x['FID']}_{x['IID']}", axis=1)
    reader_samples = set(reader.samples)
    pdata = pdata[pdata["sample"].apply(lambda x: x in reader_samples)]
    sample_order = list(pdata["sample"])
    pdata = pdata[["phenotype", "sample"] + covarcols]
    logging.info(f"Left with {pdata.shape[0]} samples after matching with VCF")

    # Get data to condition on
    if args.condition:
        logging.info(f"Loading condition on position {args.condition}")
        cond_gt = load_condition(args.vcf, args.condition, sample_order)
        pdata["condition"] = cond_gt[0]
        covarcols.append("condition")

    # Prepare output file
    outf = sys.stdout if args.out == "stdout" else open(args.out, "w")
    print_header(outf, comment_lines=[" ".join(sys.argv)])

    # Perform association tests
    logging.info(f"Performing associations with covariates: {covarcols}")
    for record in (reader.fetch(args.region) if args.region else reader):
        aaf = sum(record.aaf)
        aaf = min([aaf, 1 - aaf])
        if aaf < args.minmaf and args.minmaf != 1:
            continue

        # Analyze only STRs by default
        if len(record.REF) < MIN_STR_LENGTH:
            continue

        logging.info("Loading genotypes...")
        gts, exclude_samples = load_gt(record, sample_order, is_str=True, rmrare=args.remove_rare_str_alleles)
        pdata["GT"] = gts
        pdata["intercept"] = 1

        logging.info("Performing association...")
        assoc = perform_association(pdata, covarcols, maf=aaf, exclude_samples=exclude_samples, maxiter=args.max_iter)

        logging.info("Outputting association...")
        output_assoc(record.CHROM, record.POS, assoc, outf, assoc_type="STR")

        # Allele-based tests
        if args.allele_tests:
            logging.info("Performing allele-based tests...")
            alleles = [record.REF] + record.ALT
            for i in range(len(record.ALT) + 1):
                gts, exclude_samples = load_gt(record, sample_order, is_str=True, use_alt_num=i)
                pdata["GT"] = gts
                assoc = perform_association(pdata, covarcols, maf=aaf, exclude_samples=exclude_samples, maxiter=args.max_iter)
                output_assoc(record.CHROM, record.POS, assoc, outf, assoc_type=f"STR-alt-{alleles[i]}")

        # Allele-length-based tests (default)
        logging.info("Performing allele-length-based tests...")
        for length in set([len(record.REF)] + [len(alt) for alt in record.ALT]):
            gts, exclude_samples = load_gt(record, sample_order, is_str=True, use_alt_length=length)
            pdata["GT"] = gts
            assoc = perform_association(pdata, covarcols, maf=aaf, exclude_samples=exclude_samples, maxiter=args.max_iter)
            output_assoc(record.CHROM, record.POS, assoc, outf, assoc_type=f"STR-length-{length}")

    logging.info("Association testing complete.")

if __name__ == "__main__":
    main()
