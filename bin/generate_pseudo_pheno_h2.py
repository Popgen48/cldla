import sys
from random import shuffle

def generate_pseudo_pheno(phe_file, num_dataset):
    lc_l = [] #store the line content of the file excluding the col_pheno
    pheno_l = [] #store the original file 
    with open(phe_file) as source:
        for line in source:
            line = line.rstrip().split()
            pheno_l.append(line[2])
            del line[2]
            lc_l.append(line)
    for n in range(int(num_dataset)):
        shuffle(pheno_l)
        with open(f"dataset_{n}.phe","w") as dest:
            for i,v in enumerate(lc_l):
                v.insert(2,pheno_l[i])
                dest.write(" ".join(v)+"\n")
                del v[2] # avoid the aliasing issue

if __name__ == "__main__":
    generate_pseudo_pheno(sys.argv[1], sys.argv[2])
