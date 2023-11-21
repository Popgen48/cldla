import pysam
import sys
import random

vcf_file_path = sys.argv[1]
chrom = sys.argv[2]
window_size = int(sys.argv[3])
n_random_window = int(sys.argv[4])


vcf_file = pysam.VariantFile(vcf_file_path, 'r')

window_count = sum(1 for _ in vcf_file) - window_size + 1

vcf_file.close()

output_file_path = f'{chrom}_window_counts.txt'
random_window_file = f'{chrom}_random_windows_{n_random_window}.txt'


with open(output_file_path, 'w') as output_file:
    outprefix = 2
    for n_window in range(1,window_count+1):
        if n_window == 1:
            while outprefix <= (window_size/2+1):
                output_file.write(f'{chrom} {n_window} {outprefix}\n')
                outprefix += 1
        else:
            output_file.write(f'{chrom} {n_window} {outprefix}\n')
            outprefix += 1

n_random_l = random.sample(list(range(window_count)),n_random_window)

with open(random_window_file,'w') as output_file:
    for i,v in enumerate(n_random_l):
        output_file.write(f"{chrom}.{v}\n")

print(f'window counts for {chrom} saved to {output_file_path}')
print(f'random windows for {chrom} saved to {random_window_file}')
