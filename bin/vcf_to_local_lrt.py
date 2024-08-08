import sys
import os
import os.path
import distutils
import re
import pysam
import subprocess
import argparse
import statistics
import random
from multiprocessing import Pool
from filter_vcf import make_sample_list
from random import shuffle


abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
asreml_path = distutils.spawn.find_executable("asreml")
cldla_snp = distutils.spawn.find_executable("cldla_snp")
cldla_snp = cldla_snp if cldla_snp else f'{dname}/cldla_snp'
bend = distutils.spawn.find_executable("bend")
bend = bend if bend else f'{dname}/bend'
ginverse= distutils.spawn.find_executable("ginverse")
ginverse = ginverse if ginverse else f'{dname}/ginverse'
blupf90= distutils.spawn.find_executable("blupf90+")
blupf90 = blupf90 if blupf90 else f'{dname}/blupf90+'


g_list_geno_dict = [{}]
g_list_homo_dict = [{}]
g_list_pos_dict = [{}]
g_job_list = []


class AsremlMethods:

    def __init__(self):
        self.rm_file_pattern = ".{ask,msv,res,rsv,tmp,tsv,veo,vvp,yht}"

    def prepare_params(self, infile, numdiplo, grm, prefix, model):
        outfile_list = [f"{prefix}.as"] if model == "h1" else [grm.rstrip(".giv") + ".as"]
        for i, v in enumerate(outfile_list):
            ilc = 0
            with open(v, "w") as dest:
                with open(infile) as source:
                    for line in source:
                        line = line.rstrip()
                        if not line.startswith("!") and ilc == 0:
                            dest.write(line)
                            dest.write("\n")
                        elif ilc == 0:
                            if model == "h1":
                                dest.write(f"iDip {numdiplo} {line}")
                                dest.write("\n")
                            ilc += 1
                        elif ilc == 1:
                            dest.write(f"{grm} {line}")
                            dest.write("\n")
                            ilc += 1
                        elif ilc == 2:
                            if model == "h1":
                                dest.write(f"{prefix}.giv {line}")
                                dest.write("\n")
                            ilc += 1
                        elif ilc == 3:
                            if model == "h1" and i == 0:
                                dest.write(f"{prefix}.dat {line}")
                            elif model == "h1" and i == 1:
                                dest.write(f"{prefix}.perm.dat {line}")
                            else:
                                dest.write(f"{prefix}.h0.dat {line}")
                            dest.write("\n")
                            ilc += 1
                        elif ilc == 4:
                            if model == "h1":
                                if "iDip" not in line:
                                    print("iDip should be present in the parameter file")
                                    sys.exit(1)
                            else:
                                line = line.replace(" giv2(iDip) !r", "")
                                line = line.replace(" giv(iDip,2) !r", "")
                            dest.write(f"{line}")
                            dest.write("\n")

    def prepare_new_dat_file(self, phe):
        with open(f"{phe}.h0.dat", "w") as dest:
            with open(phe) as source:
                header = True
                for line in source:
                    if header:
                        header = False
                    else:
                        dest.write(line)

    def run_h0(self, params_file, diplo, grm, phe, model):
        prefix = grm.rstrip(".giv")
        self.prepare_params(params_file, diplo, grm, phe, model)
        self.prepare_new_dat_file(phe)
        command = f"{asreml_path} -NS5 {prefix}.as && rm {prefix}{self.rm_file_pattern}"
        subprocess.call([command], shell=True, executable='/bin/bash')
        log_l = self.extract_logl(f"{prefix}.asr")
        return log_l

    def run_asreml(self, list_i):
        u = util()
        prefix, h0_logl, mid_win_point, grm = list_i
        command = f"{asreml_path} -NS5 {prefix}.as && rm {prefix}{self.rm_file_pattern}"
        subprocess.call([command], shell=True, executable='/bin/bash')
        h1_logl = self.extract_logl(f"{prefix}.asr")
        if h0_logl != "na" and h1_logl != "na":
            return prefix, mid_win_point, -2 * (h0_logl - h1_logl)

    def extract_logl(self, asr):
        conv_logl = "na"
        last_line = "na"
        with open(asr) as source:
            for line in source:
                line = line.rstrip()
                if "LogL" in line:
                    pattern = re.compile(r"LogL=(\s){0,5}([0-9\-\.]+)")
                    match = re.findall(pattern, line)
                    if len(match) > 0:
                        conv_logl = round(float(match[0][1]), 4)
                last_line = line
        if last_line == "na":
            print(f"last line of {asr} could not be read")
        elif "LogL Converged" in last_line:
            pass
        else:
            conv_logl = "na"
        return conv_logl

    def generate_permutation_param(self, prefix):
        with open(f"{prefix}.h0.perm.as", "w") as dest1:
            with open(f"{prefix}.h1.perm.as", "w") as dest2:
                with open(f"{prefix}.as") as source:
                    for line in source:
                        line = line.rstrip()
                        if f"{prefix}.dat" in line:
                            dest1.write(line.replace(f"{prefix}.dat", f"{prefix}.h0.perm.dat") + "\n")
                            dest2.write(line.replace(f"{prefix}.dat", f"{prefix}.h1.perm.dat") + "\n")
                        elif line.startswith("iDip"):
                            dest2.write(line + "\n")
                        elif "giv" in line and "iDip" in line:
                            dest2.write(line + "\n")
                            line = line.replace(" giv2(iDip) !r", "")
                            line = line.replace(" giv(iDip,2) !r", "")
                            dest1.write(line + "\n")
                        else:
                            dest1.write(line + "\n")
                            dest2.write(line + "\n")

    def run_asreml_permutation(self, prefix):
        suffix_dict = {"h0": f"{prefix}.h0.perm.as", "h1": f"{prefix}.h1.perm.as"}
        log_l_dict = {}
        for k, v in suffix_dict.items():
            command = f"{asreml_path} -NS5 {v} && rm {prefix}.{k}.perm{self.rm_file_pattern}"
            subprocess.call([command], shell=True, executable='/bin/bash')
            log_l_dict[k] = self.extract_logl(f"{prefix}.{k}.perm.asr")
        return -2 * (log_l_dict["h0"] - log_l_dict["h1"]) if "na" not in log_l_dict.values() else None

    def process_permutation_window(self, input_list):
        u = util()
        prefix, midpoint, grm = input_list
        u.generate_pseudo_pheno(prefix, 1)
        self.generate_permutation_param(prefix)
        lrt_value = self.run_asreml_permutation(prefix)
        return prefix, midpoint, lrt_value


class Blupf90Methods:

    def prepare_params(self, infile, numdiplo, grm, prefix, model):
        outfile_list = [prefix + ".params"]
        for i, v in enumerate(outfile_list):
            with open(v, "w") as dest:
                if i == 0:
                    dest.write(f"DATAFILE\n{prefix}.dat\n")
                else:
                    dest.write(f"DATAFILE\n{prefix}.perm.dat\n")
                with open(infile) as source:
                    rewrite_line = False
                    total_effect = 0
                    last_effect = None
                    param_dict = {
                        "NUMBER_OF_EFFECTS": False,
                        "OBSERVATION(S)": False,
                        "EFFECTS:": False,
                    }
                    for line in source:
                        line = line.rstrip()
                        if (not line in param_dict and not rewrite_line) or len(line) == 0:
                            dest.write(line)
                            dest.write("\n")
                        elif rewrite_line:
                            if line != "RANDOM_RESIDUAL VALUES":
                                pattern = re.compile(r"([0-9]+)")
                                num_param = int(re.findall(pattern, line)[0])
                                if not param_dict["EFFECTS:"]:
                                    if last_param == "NUMBER_OF_EFFECTS":
                                        if model == "h0":
                                            num_param = num_param
                                        else:
                                            num_param = num_param + 1
                                    else:
                                        num_param = num_param - 1
                                    dest.write(str(num_param) + "\n")
                                    rewrite_line = False
                                    param_dict[last_param] = False
                                else:
                                    num_param = num_param if total_effect == 0 else num_param - 1
                                    dest.write(str(num_param) + " " + " ".join(line.split()[1:]) + "\n")
                                    total_effect += 1
                                    last_effect = num_param
                            else:
                                if model == "h1":
                                    dest.write(f"{last_effect+2} {numdiplo} cross")
                                    dest.write("\n")
                                dest.write(line + "\n")
                                rewrite_line = False
                        else:
                            param_dict[line] = True
                            rewrite_line = True
                            last_param = line
                            dest.write(line + "\n")
                dest.write(f"RANDOM_GROUP\n1\nRANDOM_TYPE\nuser_file\nFILE\n{grm}\n(CO)VARIANCES\n1.0\n")
                if model == "h1":
                    dest.write(
                        f"RANDOM_GROUP\n{total_effect+1}\nRANDOM_TYPE\nuser_file\nFILE\n{prefix}.giv\n(CO)VARIANCES\n0.5"
                    )
                    dest.write("\n")
                dest.write(f"OPTION method VCE")
                dest.write("\n")
                dest.write(f"OPTION maxrounds 100")
                dest.write("\n")

    def extract_logl(self, log_file, ori_window):
        pattern = re.compile(r"-2logL([^0-9]+)([0-9.]+)")
        log_l = "na"
        g_negative = False
        with open(log_file) as source:
            for line in source:
                line = line.rstrip()
                if "-2logL" in line:
                    match = re.findall(pattern, line)
                    log_l = float(match[0][1])
                if "G not positive definite" in line:
                    g_negative = True
                if not ori_window:
                    if "delta convergence=" in line:
                        if "0.000000000" in line:
                            g_negative = True
        if g_negative:
            log_l = "na"
        return log_l

    def prepare_datfile(self, phe, prefix):
        header = True
        with open(prefix + ".dat", "w") as dest:
            with open(phe) as source:
                for line in source:
                    if header == True:
                        header = False
                    else:
                        line = line.rstrip().split()
                        del line[1]
                        dest.write(" ".join(line) + "\n")

    def run_h0(self, params_file, diplo, grm, phe, model):
        prefix = grm.rstrip(".giv")
        self.prepare_params(params_file, diplo, grm, prefix, model)
        self.prepare_datfile(phe, prefix)
        command = f"mkdir {prefix} && cp {prefix}.dat {prefix}.params {grm} {prefix}/ && cd {prefix} && {blupf90} {prefix}.params >& ../{prefix}.log && cd .. && rm -r {prefix}"
        subprocess.call([command], shell=True,executable='/bin/bash')
        log_l = self.extract_logl(f"{prefix}.log", True)
        return log_l

    def create_permutation_params(self, prefix, grm_prefix):
        with open(f"{prefix}.h0.perm.params", "w") as dest_h0:
            with open(f"{prefix}.h1.perm.params", "w") as dest_h1:
                replace_line = False
                with open(f"{grm_prefix}.params") as source:
                    for line in source:
                        line = line.rstrip()
                        if line.startswith("DATAFILE"):
                            dest_h0.write(f"{line}\n")
                            replace_line = True
                        elif replace_line:
                            line = line.replace(f"{grm_prefix}.dat", f"{prefix}.h0.perm.dat")
                            dest_h0.write(f"{line}\n")
                            replace_line = False
                        else:
                            dest_h0.write(f"{line}\n")
                with open(f"{prefix}.params") as source:
                    for line in source:
                        line = line.rstrip()
                        if line.startswith("DATAFILE"):
                            dest_h1.write(f"{line}\n")
                            replace_line = True
                        elif replace_line:
                            line = line.replace(f"{prefix}.dat", f"{prefix}.h1.perm.dat")
                            dest_h1.write(f"{line}\n")
                            replace_line = False
                        else:
                            dest_h1.write(f"{line}\n")

    def run_blupf90(self, list_i):
        u = util()
        prefix, h0_logl, mid_win_point, grm = list_i
        win_ginv = prefix.rstrip(".perm") + ".giv"
        command = f"mkdir {prefix} && cp {prefix}.{{dat,params}} ./{prefix} && cp {win_ginv} {prefix}/ && cp {grm} {prefix}/ && cd {prefix} && {blupf90} {prefix}.params >& ../{prefix}.log && cd .. && rm -r {prefix}"
        subprocess.call([command], shell=True, executable='/bin/bash')
        h1_logl = self.extract_logl(f"{prefix}.log", True)
        if h0_logl != "na" and h1_logl != "na":
            return prefix, mid_win_point, h0_logl - h1_logl

    def vce_permutation_h0(self, prefix, grm):
        command = f"mkdir {prefix}_h0_perm && cp {prefix}.h0.perm.{{dat,params}} ./{prefix}_h0_perm/ && cp {prefix}.giv {prefix}_h0_perm/ && cp {grm} {prefix}_h0_perm/ && cd {prefix}_h0_perm && {blupf90} {prefix}.h0.perm.params >& ../{prefix}.h0.perm.log && cd .. && rm -r {prefix}_h0_perm"
        subprocess.call([command], shell=True,executable='/bin/bash')
        h0_logl = self.extract_logl(f"{prefix}.h0.perm.log", False)
        return h0_logl

    def vce_permutation_h1(self, prefix, grm):
        command = f"mkdir {prefix}_h1_perm && cp {prefix}.h1.perm.{{dat,params}} ./{prefix}_h1_perm/ && cp {prefix}.giv {prefix}_h1_perm/ && cp {grm} {prefix}_h1_perm/ && cd {prefix}_h1_perm && {blupf90} {prefix}.h1.perm.params >& ../{prefix}.h1.perm.log && cd .. && rm -r {prefix}_h1_perm"
        subprocess.call([command], shell=True,executable='/bin/bash')
        h1_logl = self.extract_logl(f"{prefix}.h1.perm.log", False)
        return h1_logl

    def process_permutation_window(self, list_i):
        u = util()
        prefix, mid_win_point, grm = list_i
        grm_prefix = grm.rstrip(".giv")  # replace the data file in H0 with grm_prefix
        u.generate_pseudo_pheno(prefix, 1)  # this will create h0 and h1 dataset file
        self.create_permutation_params(prefix, grm_prefix)
        h0_logl = self.vce_permutation_h0(prefix, grm)
        h1_logl = self.vce_permutation_h1(prefix, grm)
        if h0_logl != "na" and h1_logl != "na":
            return prefix, mid_win_point, h0_logl - h1_logl


class util:

    def prepare_dat_file(self, hap_path, pheno_file, prefix, tool):
        indi_diplo_dict = {}
        pheno_col = 0
        with open(prefix + ".dat", "w") as dest:
            with open(hap_path) as source:
                lc = 0
                for line in source:
                    diplo = line.rstrip().split()[3]
                    indi_diplo_dict[lc] = diplo
                    lc += 1
            with open(pheno_file) as source:
                lc = 0
                header = True
                for line in source:
                    if header:
                        header = False
                    else:
                        line = line.rstrip().split()
                        pheno_col = len(line)
                        if tool == "blupf90":
                            del line[1]
                        line.append(indi_diplo_dict[lc])
                        dest.write(" ".join(line))
                        dest.write("\n")
                        lc += 1
        if pheno_col == 0:
            print("ERROR: phenotype column not set")
            sys.exit(1)
        return pheno_col

    def get_HAP(self, hap_path, sample_genotypes):
        samples = {}
        string_ids = {}  # dictionary to store counts
        id1 = 1  # for individual haplotypes
        id2 = 1  # for combined haplotypes i.e. diplotype
        max_d = 0

        with open(hap_path, "w") as file:
            i = 1
            for key, value in sample_genotypes.items():
                str1 = []
                str2 = []
                for v in value:
                    str1.append(str(v[0] + 1))
                    str2.append(str(v[1] + 1))

                str1 = "".join(str1)
                str2 = "".join(str2)

                samples[key] = [str1, str2]

                combined_substr = str1 + str2
                combined_sbstr_rev = str2 + str1

                if str1 not in string_ids:
                    string_ids[str1] = id1
                    id1 += 1
                h1 = string_ids[str1]
                if str2 not in string_ids:
                    string_ids[str2] = id1
                    id1 += 1
                h2 = string_ids[str2]
                if combined_substr not in string_ids:
                    if combined_sbstr_rev in string_ids:
                        combined_substr = combined_sbstr_rev
                    else:
                        string_ids[combined_substr] = id2
                        id2 += 1
                d = string_ids[combined_substr]
                max_d = d if d > max_d else max_d
                str1 = " ".join(list(str1))
                str2 = " ".join(list(str2))
                file.write(f"{i} {h1} {h2} {d} {str1} {str2}\n")  # tab delimited for readability
                i += 1
        return len(samples), max_d

    def get_MAP(self, map_path, positions, hzgys):
        with open(map_path, "w") as file:
            for key in hzgys.keys():
                file.write(f"{key} {positions[key]} {hzgys[key]}\n")

    def get_PAR(self, par_path, window_size, window_number, n_samples):
        with open(par_path, "w") as file:
            file.write(
                f"100\n100\n{window_size+1}\n{window_number if window_number <= int(window_size/2) else int(window_size/2)+1}\n{n_samples}"
            )

    def get_homozygosity(self, list):
        return (list[0] ** 2 + list[1] ** 2) / (sum(list) ** 2)

    def get_ac(self, record, sample_list):
        ac = [0, 0]
        for val in sample_list:
            if record.samples[val]["GT"][0] != None:
                ac[record.samples[val]["GT"][0]] += 1
            if record.samples[val]["GT"][1] != None:
                ac[record.samples[val]["GT"][1]] += 1
        return ac

    def generate_pseudo_pheno(self, phe_prefix, num_dataset):
        lc_l = []  # store the line content of the file excluding the col_pheno
        pheno_l = []  # store the phenotype of original file
        outfile = phe_prefix + ".h1.perm.dat"
        infile = phe_prefix + ".dat"
        h0_pheno = phe_prefix + ".h0.perm.dat"
        with open(infile) as source:
            for line in source:
                line = line.rstrip().split()
                pheno_l.append(line[-2])
                del line[-2]
                lc_l.append(line)
        for n in range(int(num_dataset)):
            shuffle(pheno_l)
            with open(h0_pheno, "w") as dest_h0:
                with open(outfile, "w") as dest:
                    for i, v in enumerate(lc_l):
                        v.insert(-1, pheno_l[i])
                        dest.write(" ".join(v) + "\n")
                        dest_h0.write(" ".join(v[:-1]) + "\n")
                        # del v[-2]  # avoid the aliasing issue


class VcfToLrt:

    def read_vcf(
        self,
        vcf_path,
        chromosome,
        window_size,
        num_cores,
        grm,
        pheno_file,
        param_file,
        tool,
        num_perm,
        store,
        output_dir,
        outprefix,
    ):
        num_perm = int(num_perm)
        self.u = util()  # object containing general method
        self.asr = AsremlMethods()  # object containing asreml method
        self.blp = Blupf90Methods()  # object contaning blupf90+ method

        window_size = int(window_size)  # window size in terms of number of SNPs
        vcf = pysam.VariantFile(vcf_path)  # read vcf file
        sample_list = make_sample_list(
            pheno_file
        )  # data file of phenotypes-->sample_index, sample_id, effect1 (x1), effect2 (x2),observation (y); y should always be the last column
        dataset = outprefix
        window_number = 1
        window_process_list = []  # list will collect the parameters to run asreml and blupf90
        window_perm_process_list = (
            []
        )  # list will collect the parameters to run asreml and blupf90 on permutated dataset
        store_list = []  # store the csv strings for re-analysis
        # following condition will calculate the log-likelihood value of the null hypothesis for real data
        if tool == "asreml":
            h0_logl = self.asr.run_h0(param_file, "na", grm, pheno_file, "h0")
        else:
            h0_logl = self.blp.run_h0(param_file, "na", grm, pheno_file, "h0")
        # copy ginv of grm file to the output directory location and write the output file at the same location
        if store:
            # shutil.copy(f"{grm}", f"{output_dir}/")
            prefix = grm.rstrip(".giv")
            param_ext = "as" if tool == "asreml" else "params"
            store_list.append(
                f"{prefix},{output_dir}/{grm},{output_dir}/{pheno_file},{output_dir}/{prefix}.{param_ext}"
            )  # write the window record for H0 hypothesis
            cp_command = f"cp {prefix}.{{giv,{param_ext}}} {pheno_file} {output_dir}/ "
            subprocess.call([cp_command], shell=True,executable='/bin/bash')
        # Iterate through VCF records
        for i, record in enumerate(vcf):
            last_ele = -1
            is_window_size = False if len(list(g_list_homo_dict[last_ele].keys())) < window_size else True
            ac = self.u.get_ac(record, sample_list)  # calculate tuple of allele count for each position
            homozygosity = self.u.get_homozygosity(ac)  # calculate homozygosity
            if not record.id:
                record.id = f"{record.chrom}_{record.pos}"
            # the logic between line number 385 and 401 is this:
            # why? creating overlapping window with sliding window should not read the same record twice
            # example window size = 2
            # g_list_homo_dict = [{"snp1":0.12,"snp2":0.15},{"snp1":0.12}]
            # g_list_pos_dict = [{"snp1":1.12345,"snp2":1.67892},{"snp1":1.12345}]
            # g_list_geno_dict = [{"sample1":[(0,1),(1,1)],"sample2":[(1,0),(0,0)]},{"sample1":[(0,1)],"sample2":[(1,0)]}]
            # keep adding SNP's position, homozygosity and genotypes until it reaches the first record or reaches the window with size==2
            while not is_window_size and not len(g_list_homo_dict) + last_ele == 0:
                g_list_homo_dict[last_ele][record.id] = homozygosity
                g_list_pos_dict[last_ele][record.id] = record.pos / 1e6

                # Iterate through samples
                for sample in sample_list:
                    sample_values = record.samples[sample]["GT"]
                    if sample not in g_list_geno_dict[last_ele]:
                        g_list_geno_dict[last_ele][sample] = []
                    g_list_geno_dict[last_ele][sample].append(sample_values)

                last_ele += -1

                if len(g_list_homo_dict) > 1:
                    is_window_size = False if len(list(g_list_homo_dict[last_ele].keys())) < window_size else True
                else:
                    is_window_size = True
            g_list_homo_dict.append({record.id: homozygosity})
            g_list_pos_dict.append({record.id: record.pos / 1e6})
            g_list_geno_dict.append({})
            for sample in sample_list:
                sample_values = record.samples[sample]["GT"]
                g_list_geno_dict[-1][sample] = [sample_values]
            # it is always the first element that should reach the user-defined window size, see explanation from line 378--384
            if len(g_list_homo_dict[1]) == window_size:
                sample_genotypes = g_list_geno_dict[1]
                positions = g_list_pos_dict[1]
                hzgys = g_list_homo_dict[1]
                window_number += 1
                if window_number == 2:
                    for i in range(2, int(window_size / 2) + 2):
                        win_min_point = statistics.median(list(positions.values())[:window_number])
                        window_process_list.append(
                            (
                                f"{dataset}.{chromosome}.{window_number}",
                                float(h0_logl),
                                win_min_point,
                                None if tool == "asreml" else grm,
                            )
                        )
                        g_job_list.append(
                            (
                                dataset,
                                chromosome,
                                sample_genotypes,
                                positions,
                                hzgys,
                                window_number,
                                window_size,
                                pheno_file,
                                param_file,
                                tool,
                                grm,
                            )
                        )
                        window_number += 1
                    window_number -= 1
                else:
                    win_min_point = statistics.median(list(positions.values()))
                    window_process_list.append(
                        (
                            f"{dataset}.{chromosome}.{window_number}",
                            float(h0_logl),
                            win_min_point,
                            None if tool == "asreml" else grm,
                        )
                    )
                    g_job_list.append(
                        (
                            dataset,
                            chromosome,
                            sample_genotypes,
                            positions,
                            hzgys,
                            window_number,
                            window_size,
                            pheno_file,
                            param_file,
                            tool,
                            grm,
                        )
                    )
                if num_perm > 0:
                    window_perm_process_list.append(
                        (
                            f"{dataset}.{chromosome}.{window_number}",
                            win_min_point,
                            None if tool == "asreml" else grm,
                        )
                    )
                del g_list_geno_dict[1]
                del g_list_pos_dict[1]
                del g_list_homo_dict[1]
                if len(g_job_list) == int(num_cores):
                    with Pool(processes=len(g_job_list)) as pool:
                        pool.map(self.create_ginverse, g_job_list, 1)
                    if tool != "asreml":
                        with Pool(processes=len(window_process_list)) as pool:
                            results = pool.map(self.blp.run_blupf90, window_process_list, 1)
                        with open(f"{outprefix}_{chromosome}_results.txt", "a") as dest:
                            for result in results:
                                if result:
                                    dest.write(f"{chromosome} {result[0]} {result[1]} {result[2]}\n")
                        del window_process_list[:]
                    for i, v in enumerate(g_job_list):
                        l_prefix = f"{v[0]}.{v[1]}.{v[5]}"
                        store_list.append(l_prefix)
                    del g_job_list[:]
        if len(g_job_list) > 0:
            with Pool(processes=len(g_job_list)) as pool:
                pool.map(self.create_ginverse, g_job_list, 1)
        if len(window_process_list) > 0:
            if tool == "asreml":
                with Pool(processes=1) as pool:
                    results = pool.map(self.asr.run_asreml, window_process_list, 1)
                with open(f"{outprefix}_{chromosome}_results.txt", "a") as dest:
                    for result in results:
                        if result:
                            dest.write(f"{chromosome} {result[0]} {result[1]} {result[2]}\n")
            else:
                with Pool(processes=len(window_process_list)) as pool:
                    results = pool.map(self.blp.run_blupf90, window_process_list, 1)
                with open(f"{outprefix}_{chromosome}_results.txt", "a") as dest:
                    for result in results:
                        if result:
                            dest.write(f"{chromosome} {result[0]} {result[1]} {result[2]}\n")
            for i, v in enumerate(g_job_list):
                l_prefix = f"{v[0]}.{v[1]}.{v[5]}"
                store_list.append(l_prefix)
        if num_perm > 0:
            window_perm_process_list = random.sample(window_perm_process_list, int(num_perm))
            if tool == "asreml":
                with Pool(processes=1) as pool:
                    results_perm = pool.map(self.asr.process_permutation_window, window_perm_process_list, 1)
            else:
                with Pool(processes=int(num_cores)) as pool:
                    results_perm = pool.map(self.blp.process_permutation_window, window_perm_process_list, 1)
            with open(f"{outprefix}_{chromosome}_perm_results.txt", "a") as dest:
                for result in results_perm:
                    if result:
                        dest.write(f"{chromosome} {result[0]} {result[1]} {result[2]}\n")
        if store:
            dest = open(f"{outprefix}_{chromosome}_window_info.csv", "w")
        param_ext = "as" if tool == "asreml" else "params"
        for i, v in enumerate(store_list):
            if store:
                if i == 0:
                    dest.write(f"{v}\n")
                else:
                    dest.write(f"{v},{output_dir}/{v}.giv,{output_dir}/{v}.dat,{output_dir}/{v}.{param_ext}\n")
                    cp_command = f"cp {v}.{{giv,{param_ext},dat}} {output_dir}/"
                    subprocess.call([cp_command], shell=True,executable='/bin/bash')
                    rm_command = f"rm {v}.{{giv,{param_ext},dat}}"
                    subprocess.call([rm_command], shell=True,executable='/bin/bash')
            else:
                rm_command = f"rm {v}.{{giv,{param_ext},dat}}"
                subprocess.call([rm_command], shell=True,executable='/bin/bash')
        rm_command = f"rm *.perm.*"
        subprocess.call([rm_command], shell=True,executable='/bin/bash')
        vcf.close()

    def create_ginverse(self, input_list):

        (
            dataset,
            chromosome,
            sample_genotypes,
            positions,
            hzgys,
            window_number,
            window_size,
            pheno_file,
            params_file,
            tool,
            grm,
        ) = input_list
        prefix = f"{dataset}.{chromosome}.{window_number}"
        hap_path = f"{prefix}.Hap"
        map_path = f"{prefix}.Map"
        par_path = f"{prefix}.par"

        n_samples, max_d = self.u.get_HAP(hap_path, sample_genotypes)
        self.u.get_MAP(map_path, positions, hzgys)
        self.u.get_PAR(par_path, window_size, window_number, n_samples)

        pheno_col = self.u.prepare_dat_file(hap_path, pheno_file, prefix, tool)

        if tool == "asreml":
            self.asr.prepare_params(params_file, max_d, grm, prefix, "h1")  # prepare parameter file for asreml
        else:
            self.blp.prepare_params(params_file, max_d, grm, prefix, "h1")  # prepare parameter file for blupf90

        command = (
            f"{cldla_snp} {prefix} && {bend} {prefix}.grm {prefix}.B.grm && {ginverse} {max_d} {prefix}.B.grm {prefix}.giv"
            #+ "&& rm "
            #+ prefix
            #+ ".{Hap,Map,par,grm,B.grm}"
        )

        subprocess.call([command], shell=True, executable='/bin/bash')

        print(f"Generated .dat and .giv for {window_number}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="python script to run cLDLA on a single chromosome")

    parser.add_argument("-v", "--vcf", metavar="String", help="input phased vcf", required=True)
    parser.add_argument(
        "-r",
        "--chr",
        metavar="String",
        help="chromosome id as mentioned in vcf",
        required=True,
    )
    parser.add_argument(
        "-w",
        "--window_size",
        metavar="Int",
        help="window size",
        default=40,
        required=False,
    )
    parser.add_argument(
        "-c",
        "--num_cpus",
        metavar="Int",
        help="number of windows to be run in parallel",
        default=8,
        required=False,
    )
    parser.add_argument(
        "-g",
        "--grm",
        metavar="String",
        help="relationship matrix based on entire chromosome",
        required=True,
    )
    parser.add_argument("-p", "--pheno", metavar="String", help="phenotype file", required=True)
    parser.add_argument("-a", "--params", metavar="String", help="parameter file", required=True)
    parser.add_argument(
        "-t",
        "--tool",
        default="asreml",
        const="asreml",
        nargs="?",
        choices=["asreml", "blupf90"],
        help="tools for variance component estimation --> asreml or blupf90 (default: %(default)s)",
    )

    parser.add_argument("-o", "--outprefix", help="output prefix", default="cldla", required=False)
    parser.add_argument(
        "-n",
        "--num_perm",
        metavar="Int",
        help='carry out "n" permutation to access the level of significance (True or False)',
        default=100,
        required=False,
    )
    parser.add_argument(
        "-s",
        "--store",
        help="whether or not to store ginv, grm and phenotype file for each window",
        action="store_true",
    )
    parser.add_argument("-O", "--output_dir", help="path to the output directory", default="cldla", required=False)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif args.store and not args.output_dir:
        print("when setting the tag --store the output_dir must be defined")
        sys.exit(1)
    else:
        calc_lrt = VcfToLrt()
        calc_lrt.read_vcf(
            args.vcf,
            args.chr,
            args.window_size,
            args.num_cpus,
            args.grm,
            args.pheno,
            args.params,
            args.tool,
            args.num_perm,
            args.store,
            args.output_dir,
            args.outprefix,
        )
