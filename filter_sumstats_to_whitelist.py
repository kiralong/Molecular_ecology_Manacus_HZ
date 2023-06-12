#!/usr/bin/env python3
import sys, os, argparse, random

#
# Command line options
#
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('-s', '--sumstats', required=True, help='Populations Sumstats TSV file')
    p.add_argument('-p', '--n-populations', required=False, default=1, type=int, help='Min Number of Populations')
    p.add_argument('-n', '--number-sites', required=False, default=1000, type=int, help='Number of Sites to Keep')
    p.add_argument('-e', '--hwe', required=False, action='store_true', default=False, help='Check for Sites in HWE')
    p.add_argument('-f', '--maf', required=False, default=0.0, type=float, help='Min Allele Frequency Cutoff')
    p.add_argument('-o', '--outd', required=False, default='.', help='Output directory')
    # Check input arguments
    args = p.parse_args()
    args.outd = args.outd.rstrip('/')
    if not os.path.exists(args.sumstats):
        sys.exit(f"Error: '{args.sumstats}' not found.")
    if not os.path.exists(args.outd):
        sys.exit(f"Error: '{args.outd}' not found.")
    if args.maf > 0.5:
        sys.exit(f"Error: MAF ({args.maf}) greater than 0.5.")
    if args.number_sites <= 0:
        sys.exit(f"Error: 'number-sites' ({args.number_sites}) must be non-zero positive integer.")
    return args

# Site from the Sumstats File
class SumstatsSite:
    def __init__(self, locus_id, locus_col, chrom, bp, popid, p, hwe, private):
        self.locid = locus_id
        self.locol = locus_col
        self.chrom = chrom
        self.bp    = bp
        self.popid = popid
        self.p     = p
        self.hwe   = hwe
        self.priv  = private
    def __str__(self):
        return f'{self.locid}\t{self.locol}\t{self.chrom}\t{self.bp}\t{self.popid}\t{self.p}\t{self.hwe}\t{self.priv}'

#
# Parse Sumstats file
#
def parse_sumstats(sumstats_f, n_pops, maf, hwe):
    assert 0.0 <= maf <= 0.5
    assert type(hwe) is bool
    max_hwe_p = 0.0
    if hwe is True:
        max_hwe_p = 0.05
    # Output
    sites_dict = dict()
    prev_snp = None
    pops_in_site = list()
    # Parse the sumstats
    for i, line in enumerate(open(sumstats_f, 'r')):
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        locid = int(fields[0])
        locol = int(fields[3])
        chrom = fields[1]
        bp = int(fields[2])
        popid = fields[4]
        p = float(fields[8])
        hwe_p = float(fields[19])
        priv = int(fields[20])
        snp = f'{locid}_{locol}'
        # Filter low maf
        if p >= 0.5:
            if 1-p < maf:
                continue
        else:
            if p < maf:
                continue
        # Filter by HWE
        if hwe_p < max_hwe_p:
            continue
        site = SumstatsSite(locid, locol, chrom, bp, popid, p, hwe_p, priv)
        # For the first snp
        if prev_snp is None:
            prev_snp = snp
        # If looking at the same SNP in a different pop just add to the list
        if snp == prev_snp:
            pops_in_site.append(site)
        # When looking at a different SNP
        else:
            prev_loc = int(prev_snp.split('_')[0])
            prev_col = int(prev_snp.split('_')[1])
            if len(pops_in_site) == n_pops:
                sites_dict.setdefault(prev_loc, {})
                sites_dict[prev_loc][prev_col] = pops_in_site
            pops_in_site = list()
            prev_snp = None
            pops_in_site.append(site)
        prev_snp = snp
    prev_loc = int(prev_snp.split('_')[0])
    prev_col = int(prev_snp.split('_')[1])
    if len(pops_in_site) == n_pops:
        sites_dict.setdefault(prev_loc, {})
        sites_dict[prev_loc][prev_col] = pops_in_site
    return sites_dict

# Sample the final sites
def sample_kept_sites(sites_dict, n_sites, outd):
    if len(sites_dict) == 0:
        return
    if n_sites > len(sites_dict):
        print(f"Warning: more sites chosen ({n_sites:,}) than loci available ({len(sites_dict):,}). Exporting {len(sites_dict):,} total sites.")
        n_sites = len(sites_dict)
    kept_loci = random.sample(sites_dict.keys(), n_sites)
    out = open(f'{outd}/snp_whitelist.tsv', 'w')
    for loc in sorted(kept_loci):
        col = random.choice(list(sites_dict[loc].keys()))
        out.write(f'{loc}\t{col}\n')
    out.close()

def main():
    args = parse_args()
    # Read sumstats and filter sites
    sites_dict = parse_sumstats(args.sumstats, args.n_populations, args.maf, args.hwe)
    nlocs = len(sites_dict)
    nsnps = 0
    for loc in sites_dict:
        for site in sites_dict[loc]:
            nsnps += 1
    print(f'Total loci kept: {nlocs:,}.')
    print(f'Total SNPs kept: {nsnps:,}.')
    # Subsample the kept sites and save to file
    sample_kept_sites(sites_dict, args.number_sites, args.outd)

# Run Code
if __name__ == '__main__':
    main()
