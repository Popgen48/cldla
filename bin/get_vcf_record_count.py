import pysam
import sys
import os

vcf_file_path = sys.argv[1]
chrom = sys.argv[2]
window_size = int(sys.argv[3])

vcf_file = pysam.VariantFile(vcf_file_path, 'r')

window_count = sum(1 for _ in vcf_file) - window_size + 1

vcf_file.close()

output_file_path = chrom+'_window_counts.txt'

# check if file already exists
if not os.path.exists(output_file_path):
    with open(output_file_path, 'w') as output_file:
        for n_window in range(1,window_count+1):
            output_file.write(f"{chrom} {n_window}\n")

print(f'window count saved to {output_file_path}')
