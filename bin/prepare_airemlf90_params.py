import sys

def prepare_params(infile, numdiplo, chrom_giv, chromwind_giv, pheno_file, model):
    outfile = ".".join(pheno_file.split(".")[:-1])+".params" if model == "h1" else chrom_giv.strip(".giv")+".params"
    lct = 0
    header = 0
    trait_column = 0
    read_observation = False
    write_random_effect = True
    read_effect_par = False
    read_number_effect = False
    read_rand_res = False
    with open(outfile,"w") as dest:
        with open(infile) as source:
            for line in source:
                if header == 0:
                    dest.write(f'DATAFILE\n{pheno_file}\n')
                    header += 1
                    dest.write(line)
                elif line.startswith('NUMBER_OF_EFFECTS'):
                    read_number_effect = True
                    dest.write(line)
                elif line.startswith('OBSERVATION'):
                    read_observation = True
                    dest.write(line)
                elif read_observation:
                    trait_column = int(line.rstrip().split()[0])
                    read_observation = False
                    dest.write(line)
                elif line.startswith("EFFECTS"):
                    lct +=1
                    read_effect_par = True
                    dest.write(line)
                elif line.startswith("RANDOM_RESIDUAL"):
                    read_effect_par = False
                    read_rand_res = True
                    if numdiplo != "na":
                        diplo_idx = lct + 1 if trait_column != lct+1 else lct+2
                        dest.write(f'{diplo_idx} {numdiplo} cross\n')
                    dest.write(line)
                elif read_effect_par:
                    lct += 1
                    dest.write(line)
                elif line.startswith('OPTION'):
                    if write_random_effect:
                        dest.write(f'RANDOM_GROUP\n1\nRANDOM_TYPE\nuser_file\nFILE\n{chrom_giv}\n(CO)VARIANCES\n1.0\n')
                        if numdiplo != "na":
                            dest.write(f'RANDOM_GROUP\n{lct}\nRANDOM_TYPE\nuser_file\nFILE\n{chromwind_giv}\n(CO)VARIANCES\n0.5\n')
                        dest.write("\n")
                        write_random_effect = False
                    dest.write(line)
                elif read_number_effect:
                    read_number_effect = False
                    if model == "h0":
                        dest.write(str(int(line.rstrip().split()[0])-1)+"\n")
                    else:dest.write(line)
                elif read_rand_res:
                    dest.write(line)
                    read_rand_res = False
                else:
                    dest.write(line)

if __name__ == "__main__":
    prepare_params(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
