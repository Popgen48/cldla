import sys

def prepare_params(infile, numdiplo, chrom_giv, chromwind_giv, pheno_file):
    outfile = ".".join(chromwind_giv.split(".")[:-1])+".as"
    ilc = 0
    with open(outfile,"w") as dest:
        with open(infile) as source:
            for line in source:
                line = line.rstrip()
                if not line.startswith("!"):
                    dest.write(line)
                elif ilc == 0:
                    dest.write(f"iDip {numdiplo} {line}")
                    ilc += 1
                elif ilc ==1:
                    dest.write(f"{chrom_giv} {line}")
                    ilc += 1
                elif ilc == 2:
                    dest.write(f"{chromwind_giv} {line}")
                    ilc += 1
                elif ilc == 3:
                    dest.write(f"{pheno_file} {line}")
                    ilc += 1
                elif ilc == 4:
                    if "iDip" not in line:
                        print("iDip should be present in the parameter file")
                        sys.exit(1)
                    else:
                        dest.write(f"{line}")
                dest.write("\n")

if __name__ == "__main__":
    prepare_params(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
