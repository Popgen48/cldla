import sys
import pandas as pd
from extract_chromosome_ids import read_input_file

chrom_vcf_idx = sys.argv[1]
outprefix = sys.argv[2]
process = sys.argv[3]
n_perm = sys.argv[4]
chrom_files = sys.argv[5:]

global_list = []

chrom_replace_d, chrom_is_string, chrom_is_chr = read_input_file(chrom_vcf_idx)

for c_file in chrom_files:
    with open(c_file) as source:
        for line in source:
            line = line.rstrip()
            if not "None" in line:
                if chrom_is_string or chrom_is_chr:
                    new_record = line.replace(line.split()[0], chrom_replace_d[line.split()[0]]).split()
                else:
                    new_record = line.split()
                global_list.append(new_record)

df = pd.DataFrame(global_list, columns=["chrom", "window_name", "cord", "LRT_values"])

if int(n_perm) == 0:
    chrom_list = df["chrom"].unique().tolist()
    with open(f"{outprefix}_chrom_threshold.txt","w") as dest:
        for chrom in chrom_list:
            dest.write(f"{chrom}\t18.00\n") # ---> this is the default value for each local threshold
else:
    df = df.drop("window_name", axis=1)
    df["chrom"] = pd.to_numeric(df["chrom"])
    df["cord"] = pd.to_numeric(df["cord"]) * 1000000
    df["cord"] = df["cord"].astype(int)
    df["LRT_values"] = df["LRT_values"].astype(float)
    df.loc[df["LRT_values"] < 0, "LRT_values"] = 0

    if process != "estimate":
        df = df.sort_values(["chrom", "cord"], ascending=[True, True])
        with open(f"{outprefix}_maninp.txt", "w") as dest:
            dest.write(df.to_string(header=False, index=False))
    else:
        df = df.sort_values(["chrom", "LRT_values"], ascending=[True, False])
        df1 = df.groupby("chrom").first()
        df1 = df1.reset_index()
        del df1["cord"]
        with open(f"{outprefix}_chrom_threshold.txt", "w") as dest:
            dest.write(df1.to_string(header=False, index=False))
