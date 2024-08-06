import sys

pheno_file = sys.argv[1]
tool = sys.argv[2]
extra_params = sys.argv[3]
outprefix = sys.argv[4]

extra_params_line = ""

dest = open(outprefix+"_updated.phe","w")

if tool != "asreml" and tool != "blupf90":
    print(f"the tool name should either be asreml or blupf90, provided tool {tool}")
    sys.exit(1)

with open(pheno_file) as source:
    header = True
    no_fixed_effect = False
    for line in source:
        line = line.rstrip().split()
        if header:
            if line[:2] != ["iidx", "iid"]:
                print(
                    f"the phenotype file {pheno_file} does not contain header line, header line should contain the two columns with header iidx and iid, provided header of the first two columns {line[:2]}"
                )
                sys.exit(1)
            else:
                no_fixed_effect = True if len(line)==3 and tool == "blupf90" else False
                if no_fixed_effect:
                    line.insert(2,"fix_intrp")
                header_info = {k: [] for k in line}
                header_keys = list(header_info.keys())
            header = False
        else:
            if no_fixed_effect:
                line.insert(2,"1")
            for i, v in enumerate(line):
                if v not in header_info[header_keys[i]]:#to count the level of effects, the level should be unique
                    header_info[header_keys[i]].append(v)
        dest.write(f'{" ".join(line)}\n')
dest.close()

#decide whether the fixed effect is "covariate" or "crossclassified", simple logic is--> "." is present then "covariate" otherwise "cross"

type_effect = {}
for i,v in enumerate(header_keys):
    header_values = header_info[v]
    header_cov_info = [True if "." in v else False for v in header_values]
    type_effect[v] = "cov" if True in header_cov_info else "cross"
    


ext = "as" if tool == "asreml" else "params"
dest = open(f"{outprefix}.{ext}", "w")

if extra_params != "None":
    with open(extra_params) as source:
        for line in source:
            line = line.rstrip()
            extra_params_line += f"{line}\n"


if tool == "asreml":
    dest.write(f"This line is a headline: MLMs for the trait {header_keys[-1]}\n")
    for ke in header_keys:
        if ke in ["iidx", "iid"]:
            dest.write(f"{ke} {len(header_info[ke])} {'!I' if ke == 'iidx' else '!A'}\n")
        else:
            dest.write(f"{ke} 1\n")
    dest.write(f"!I\n!LDET\n!LDET\n")
    if extra_params_line != "":
        dest.write(f"{extra_params_line}")
    dest.write(f'{header_keys[-1]} ~ mu {" ".join([k for k in header_keys[2:-1]])} !r giv(iDip,2) !r giv(iid,1)\n')
else:
    dest.write(
        f"NUMBER_OF_TRAITS\n1\nNUMBER_OF_EFFECTS\n{len(header_keys)-2}\nOBSERVATION(S)\n{len(header_keys)}\nWEIGHT(S)\n\nEFFECTS:\n"
    )
    for i, v in enumerate(header_keys):
        if i != 1 and i != len(header_keys) - 1:
            dest.write(f'{i+1} {len(header_info[v] if type_effect[v] == "cross" else "1")} {"cross" if i==0 else type_effect[v]}\n')
    if extra_params_line != "":
        dest.write(f"{extra_params_line}")
