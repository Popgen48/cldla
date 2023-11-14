import sys
import re

def convert_ginv(log_file,g_inv, wg_inv):
    with open(wg_inv, "w") as dest:
        with open(log_file) as source:
            for line in source:
                if "determinant" in line:
                    line=line.rstrip()
                    pattern=re.compile(r'(Log determinant=)([\s]{0,1}[\-0-9\.]*)(.*)')
                    matcher=re.findall(pattern,line)
                    match=matcher[0][1].lstrip()
                    dest.write(match+" "+"dense")
                    dest.write("\n")
        with open(g_inv) as source:
            for line in source:
                line=line.rstrip()
                a=line.split()
                dest.write(a[1]+" "+a[0]+" "+a[2])
                dest.write("\n")

if __name__ == "__main__":
    convert_ginv(sys.argv[1],sys.argv[2], sys.argv[3])
