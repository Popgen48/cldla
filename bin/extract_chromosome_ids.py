import sys
import re


def read_input_file(chrom_vcf_idx):
    chrom_replace_d = {}
    chrom_count = 0
    chrom_is_num = False
    chrom_is_string = False
    chrom_is_chr = False
    line_count = 0
    pattern = re.compile(r"([chr|CHR|Chr][A-Za-z]+)([0-9]+)")
    with open(chrom_vcf_idx) as source:
        for line in source:
            line = line.rstrip().split(",")
            if not line[0].isdigit():
                match = re.findall(pattern, line[0])
                if len(match) == 0:
                    chrom_count += 1
                    chrom_is_string = True
                    chrom_replace_d[line[0]] = str(chrom_count)
                else:
                    chrom_replace_d[line[0]] = match[0][1]
                    chrom_is_chr = True
            else:
                chrom_is_num = True
            line_count += 1

    if (chrom_is_num and chrom_is_string) or (chrom_is_num and chrom_is_chr) or (chrom_is_chr and chrom_is_num):
        print("the chromosome column in the input csv file contains string as well as integer chromosome id")
        sys.exit(1)
    return chrom_replace_d, chrom_is_string, chrom_is_chr
