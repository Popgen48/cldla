# Simple script to extract the values from all the samples of a given .VCF file, assign them an id and store in a .hap file
import sys
import pysam
import argparse


def make_sample_list(incl_samples):
    sample_list = []
    with open(incl_samples) as source:
        for line in source:
            line = line.rstrip().split()
            sample_list.append(line[1])
    return sample_list


def get_maf(record, sample_list):
    ac = [0, 0]
    for val in sample_list:
        if record.samples[val]["GT"][0] != None:
            ac[record.samples[val]["GT"][0]] += 1
        if record.samples[val]["GT"][1] != None:
            ac[record.samples[val]["GT"][1]] += 1
    return ac


def create_new_header(vcf_header, sample_list):
    new_header = pysam.VariantHeader()
    ordered_sample_list = []  # created this list to make sure that the order in the header and record stays the same
    for record in vcf_header.records:
        new_header.add_record(record)
    for sample in vcf_header.samples:
        if sample in sample_list:
            ordered_sample_list.append(sample)
            new_header.add_sample(sample)
    return new_header, ordered_sample_list


def read_vcf(vcf_path, incl_samples, maf_threshold, output_file):
    sample_list = make_sample_list(incl_samples)
    vcf = pysam.VariantFile(vcf_path)
    new_header, ordered_sample_list = create_new_header(vcf.header, sample_list)
    vcf_out = pysam.VariantFile(output_file, "w", header=new_header)
    vcf.subset_samples(ordered_sample_list)

    # Iterate through VCF records
    for i, record in enumerate(vcf):
        ac = get_maf(record, sample_list)
        maf = min(ac) / sum(ac)
        # calculate maf only collect the genotypes if the maf is grt than the threshold
        if maf > float(maf_threshold):
            vcf_out.write(record)
    vcf.close()
    vcf_out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="python script to filter vcf files based on maf and samples"
    )

    parser.add_argument("-v", "--vcf", metavar="String", help="input phased vcf file", required=True)
    parser.add_argument("-s", "--sample_file", metavar="String", help="samples to be included", required=True)
    parser.add_argument(
        "-m",
        "--maf",
        metavar="Float",
        help="minimum MAF required to be included inanalysis",
        default=0.05,
        required=False,
    )
    parser.add_argument("-o", "--vcf_out", metavar="String", help="output vcf file after filtering", required=True)
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        read_vcf(
            args.vcf,
            args.sample_file,
            args.maf,
            args.vcf_out,
        )
