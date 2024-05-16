import sys

phe_file = sys.argv[1]
outprefix = sys.argv[2]

phe_list = []
qcovar_list = []
covar_list = []

qcovar_col = []
covar_col = []

row_count = 0
with open(phe_file) as source:
    for line in source:
        row_count += 1
        line = line.rstrip().split()
        phe_list.append(line[1]+" "+line[1]+" "+line[-1])
        covar_t_list = [line[1],line[1]]
        qcovar_t_list = [line[1],line[1]]
        for i,v in enumerate(line[2:-1]):
            if "." in v:
                qcovar_t_list.append(v)
                if i not in qcovar_col:
                    qcovar_col.append(i)
                if i in covar_col:
                    print(f"ERROR:in column {i+2}, row {row_count} there is a mixed of float and integer values")
                    sys.exit(1)
            else:
                covar_t_list.append(v)
                if i not in covar_col:
                    covar_col.append(i)
                if i in qcovar_col:
                    print(f"ERROR:in column {i+2}, row {row_count} there is a mixed of float and integer values")
                    sys.exit(1)
        if len(qcovar_t_list)>2:
            qcovar_list.append(" ".join(qcovar_t_list[:]))
        if len(covar_t_list)>2:
            covar_list.append(" ".join(covar_t_list[:]))
        del covar_t_list[:]
        del qcovar_t_list[:]

with open(f'{outprefix}.phe',"w") as dest:
    dest.write("\n".join(phe_list))
if len(qcovar_list)>0:
    with open(f'{outprefix}.qcovar',"w") as dest:
        dest.write("\n".join(qcovar_list))
if len(covar_list)>0:
    with open(f'{outprefix}.covar',"w") as dest:
        dest.write("\n".join(covar_list))
