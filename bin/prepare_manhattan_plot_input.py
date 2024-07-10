import sys
import pandas as pd
from extract_chromosome_ids import read_input_file

chrom_vcf_idx = sys.argv[1]
outprefix = sys.argv[2]
chrom_files = sys.argv[3:]

global_list = []

chrom_replace_d, chrom_is_string, chrom_is_chr = read_input_file(chrom_vcf_idx)

for c_file in chrom_files:
    with open(c_file) as source:
        for line in source:
            line = line.rstrip()
            if chrom_is_string or chrom_is_chr:
                new_record = line.replace(line.split()[0],chrom_replace_d[line.split()[0]]).split()
            else:
                new_record = line.split()
            global_list.append(new_record)
df = pd.DataFrame(global_list, columns = ['chrom', 'window_name','cord', 'LRT_values']) 
df = df.drop('window_name', axis=1)
df["chrom"] = pd.to_numeric(df["chrom"])
df["cord"] = pd.to_numeric(df["cord"])*1000000
df = df.sort_values(['chrom', 'cord'], ascending=[True, True])

with open(f"{outprefix}_maninp.txt","w") as dest:
    dest.write(df.to_string(header=False, index=False))
