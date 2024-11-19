# Simple script to extract the values from all the samples of a given .VCF file, assign them an id and store in a .hap file
import sys
import pysam


def read_vcf(vcf_path, snp_window_size):
    nonoverlap_window_d = {}
    vcf = pysam.VariantFile(vcf_path)

    # Iterate through VCF records
    for record in vcf:
        if record.chrom not in nonoverlap_window_d:
            nonoverlap_window_d[record.chrom] = 0
        nonoverlap_window_d[record.chrom] += 1

    vcf.close()

    with open("chrom_window_record_counts.txt", "w") as dest:
        for chrom in nonoverlap_window_d:
            dest.write(f"{chrom} {nonoverlap_window_d[chrom]//int(snp_window_size)}\n")


if __name__ == "__main__":
    read_vcf(sys.argv[1], sys.argv[2])
