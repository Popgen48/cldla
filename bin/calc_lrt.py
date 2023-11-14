import sys

# takes log-likelihood file of h0 hypothesis and log-likelihood file of h1 hypothesis as inputs and generate *.lrt.out

def calc_lrt(chrom, h0_llik, h1_llik):

    #read h0 log-likelihood file first
    with open(h0_llik) as source:
        for line in source:
            line = line.rstrip().split()
            Lo = float(line[-1]) # set L0 value
    print(Lo)
    with open(f'{chrom}.lrt.out',"w") as dest:
        with open(h1_llik) as source:
            for line in source:
                print(line)
                line = line.rstrip().split()
                lrt = -2*(Lo-float(line[-1]))
                dest.write(f'{chrom} {line[1]} {lrt}')
                dest.write("\n")

if __name__ == "__main__":
    calc_lrt(sys.argv[1], sys.argv[2], sys.argv[3])
