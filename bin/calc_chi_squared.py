import sys
from scipy.stats import chi2
from math import log10


in_df = sys.argv[1]
num_snps = sys.argv[2]
out_df = sys.argv[3]

chrom_list = []
pval_lrt_d = {}

cutoff = "na"


with open(out_df + ".chisq_2log10Pval.maninp.txt", "w") as dest1:
    with open(out_df + ".lrt.chisqPval.maninp.txt", "w") as dest:
        with open(in_df) as source:
            for line in source:
                line = line.rstrip().split()
                p_val = round(float(chi2.sf(float(line[2]), 1)),8)
                dest.write(f"{line[0]} {line[1]} {line[2]} {p_val}\n")
                dest1.write(f"{line[0]} {line[1]} {-log10(p_val)}\n")
                if f"{chi2.sf(float(line[2]),1)}" not in pval_lrt_d:
                    pval_lrt_d[p_val] = line[2]
                if line[0] not in chrom_list:
                    chrom_list.append(line[0])
pval_lrt_sd = {k: pval_lrt_d[k] for k in sorted(pval_lrt_d)}

with open(num_snps) as source:
    total_win = 0
    for line in source:
        line = line.rstrip().split()
        if line[0] in chrom_list:
            total_win += int(line[1])

multi_corr_p_val = 0.05 / total_win

for k in pval_lrt_sd:
    if k <= multi_corr_p_val:
        cutoff = pval_lrt_sd[k]
with open(f"{out_df}.lrt.cutoff", "w") as dest:
    with open(f"{out_df}.pval.cutoff", "w") as dest1:
        for chrom in chrom_list:
            dest.write(f"{chrom} {cutoff}\n")
            dest1.write(f"{chrom} {-log10(multi_corr_p_val)}\n")
