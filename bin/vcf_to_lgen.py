# Simple script to extract the values from all the samples of a given .VCF file, assign them an id and store in a .hap file
import sys
import pysam

def pheno_to_list(pheno_file):
    indi_list = []
    with open(pheno_file) as source:
        for line in source:
            line = line.rstrip().split()
            indi_list.append(line[1])
    return indi_list

def read_vcf(vcf_path, chrm, indi_list):
    # dictionary to store samples
    sample_genotypes = {}
    # dictionary to store MAF
    snpid_l = []
    # dictionary to store homozygosity
    vcf = pysam.VariantFile(vcf_path)

    # Iterate through VCF records
    for record in vcf:
        if record.chrom != chrm:
            allele_base = {
                0: record.ref[0],
                1: record.ref[0] if record.alts[0] == "." else record.alts[0],
            }
            for sample in indi_list:
                genotypes_t = record.samples[sample]["GT"]
                genotype_lt = []
                if sample not in sample_genotypes:
                    sample_genotypes[sample] = []

                genotype_lt.append(allele_base[genotypes_t[0]])
                genotype_lt.append(allele_base[genotypes_t[1]])
                sample_genotypes[sample].append(genotype_lt[:])
                del genotype_lt[:]

            if not record.id:
                record.id = f"{record.chrom}_{record.pos}"

            snpid_l.append(record.id)

    vcf.close()

    return sample_genotypes, snpid_l


def write_lgen(sample_genotypes, snpid_l, outprefix):
    with open(outprefix + ".lgen", "w") as dest:
        for sample in sample_genotypes:
            genotype_l = sample_genotypes[sample]
            for i, v in enumerate(genotype_l):
                dest.write(sample + "\t" + snpid_l[i] + "\t" + "\t".join(v) + "\n")


def vcf_to_lgen(vcf_path, chrom, pheno_file, outprefix):
    indi_list = pheno_to_list(pheno_file)
    sample_genotypes, snpid_l = read_vcf(vcf_path, chrom, indi_list)
    write_lgen(sample_genotypes, snpid_l, outprefix)


if __name__ == "__main__":
    vcf_to_lgen(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
