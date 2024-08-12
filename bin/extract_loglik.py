import sys


def read_map(map_file, window_size):
    pos_list = []
    window_count = int(map_file.split(".")[-2])
    with open(map_file) as source:
        for line in source:
            line = line.rstrip().split()
            pos_list.append(float(line[1]))
    if window_count < int(window_size) / 2 + 1:
        median = pos_list[window_count]
    elif len(pos_list) % 2 == 0:
        i_mid = int(len(pos_list) / 2) - 1
        median = (pos_list[i_mid] + pos_list[i_mid + 1]) / 2
    else:
        i_mid = int(len(pos_list) / 2)
        median = pos_list[i_mid]
    return median * 1000000


def extract_log_likelihood(esr_file, map_file, out_file, chrom, window_size, tool):
    median = read_map(map_file, window_size) if map_file != "none" else "00"
    chrom_window = (
        map_file.split(".")[-3] + "." + map_file.split(".")[-2] if map_file != "none" else esr_file.split("_")[0]
    )
    chrom = chrom_window.split(".")[0] if map_file != "none" else chrom
    print(chrom_window)
    window = chrom_window.split(".")[1] if map_file != "none" else "00"
    with open(out_file, "a+") as dest:
        if tool != "blupf90+":
            with open(esr_file) as source:
                log_l = "na"
                try:
                    for line in source:
                        if "Error" not in line:
                            if "LogL=" in line:
                                line = line.rstrip().split()
                                log_l = line[2]
                            if "Finished" in line:
                                if "LogL Converged" in line and not "not" in line:
                                    dest.write(f"{chrom} {window} {median} {log_l}\n")
                                else:
                                    print(f"echidna for {chrom} {window} either not finished or not converged")
                                    sys.exit(1)
                except:
                    print(f"LogL cannot be extracted from {esr_file} ")
        else:
            with open(esr_file) as source:
                for line in source:
                    dest.write(f"{chrom} {window} {median} {line}")


if __name__ == "__main__":
    extract_log_likelihood(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
