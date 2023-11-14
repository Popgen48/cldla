import sys

def prepare_params(infile, numdiplo, pheno_file, model):
    outfile = "wombat.par"
    with open(outfile,"w") as dest:
        with open(infile) as source:
            for line in source:
                if "DATA" in line:
                    line = line.rstrip()
                    line = line+pheno_file+"\n"
                if line.startswith("Ibd"):
                    line = line.rstrip()
                    line = line+" "+numdiplo+"\n"
                dest.write(line)

if __name__ == "__main__":
    prepare_params(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
