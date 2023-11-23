import sys

def read_map(map_file, window_size):
    pos_list = []
    window_count = int(map_file.split(".")[-2])
    with open(map_file) as source:
        for line in source:
            line = line.rstrip().split()
            pos_list.append(float(line[1]))
    if window_count<int(window_size)/2+1:
        median = pos_list[window_count]
    elif len(pos_list) % 2 == 0:
        i_mid = int(len(pos_list)/2)-1
        median = ((pos_list[i_mid] + pos_list[i_mid+1])/2)
    else:
        i_mid = int(len(pos_list)/2)
        median = pos_list[i_mid]
    return median*1000000

def extract_log_likelihood(esr_file, map_file, out_file, window_size):
    median = read_map(map_file, window_size) if map_file != "none" else "00"
    chrom_window = esr_file.split(".")[-3]+"_"+esr_file.split(".")[-2] if map_file != "none" else esr_file.split("_")[0]
    with open(out_file,"a+") as dest:
        with open(esr_file) as source:
            log_l = "na"
            try:
                for line in source:
                    if "Error" not in line:
                        if "LogL=" in line:
                            line = line.rstrip().split()
                            log_l = line[2]
                        if "Finished" in line:
                            if "LogL Converged" in line and not "not"  in line:
                                dest.write(f'{chrom_window} {median} {log_l}\n')
                            else:
                                print(f'echidna for {chrom} {window} either not finished or not converged')
                                sys.exit(1)
            except:
                print(f'LogL cannot be extracted from {esr_file} ')
                #dest.write(f'{chrom_window} {median} {log_l}\n')


if __name__ == "__main__":
    extract_log_likelihood(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
