import sys
import os
import re
import pysam
import subprocess
import argparse
import numpy as np
from multiprocessing import Pool
from filter_vcf import make_sample_list

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)

g_list_geno_dict = [{}]
g_list_homo_dict = [{}]
g_list_pos_dict = [{}]
g_job_list = []


class AsremlMethods:

    def prepare_params(self, infile, numdiplo, grm, prefix, model):
        outfile = f"{prefix}.as" if model == "h1" else grm.rstrip(".giv") + ".as"
        ilc = 0
        with open(outfile, "w") as dest:
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
                        if model == "h1":
                            dest.write(f"{prefix}.dat {line}")
                        else:
                            dest.write(f"{prefix} {line}")
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

    def run_h0(self, params_file, diplo, grm, phe, model):
        prefix = grm.rstrip(".giv")
        self.prepare_params(params_file, diplo, grm, phe, model)
        rm_file_pattern = ".{ask,aov,msv,res,rsv,sln,tmp,tsv,veo,vvp,yht}"
        command = f"asreml -NS5 {prefix}.as && rm {prefix}{rm_file_pattern}"
        subprocess.call([command], shell=True)
        log_l = self.extract_logl(f"{prefix}.asr")
        return log_l

    def run_asreml(self, list_i):
        prefix, h0_logl, mid_win_point = list_i
        rm_file_pattern = ".{ask,aov,msv,res,rsv,sln,tmp,tsv,veo,vvp,yht}"
        command = f"asreml -NS5 {prefix}.as && rm {prefix}{rm_file_pattern}"
        subprocess.call([command], shell=True)
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
                    pattern = re.compile(r"LogL=(\s){0,5}([0-9\-]+)")
                    match = re.findall(pattern, line)
                    if len(match) > 0:
                        conv_logl = float(match[0][1])
                last_line = line
        if last_line == "na":
            print(f"last line of {asr} could not be read")
        elif "LogL Converged" in last_line:
            pass
        else:
            conv_logl = "na"
        return conv_logl


class Blupf90Methods:

    def prepare_params(self, infile, numdiplo, grm, prefix, model):
        param_dict = {
            "NUMBER_OF_EFFECTS": False,
            "OBSERVATION(S)": False,
            "EFFECTS:": False,
        }
        rewrite_line = False
        total_effect = 0
        with open(prefix + ".params", "w") as dest:
            dest.write(f"DATAFILE\n{prefix}.dat\n")
            with open(infile) as source:
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
                                num_param = (
                                    num_param if total_effect == 0 else num_param - 1
                                )
                                dest.write(
                                    str(num_param)
                                    + " "
                                    + " ".join(line.split()[1:])
                                    + "\n"
                                )
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
            dest.write(
                f"RANDOM_GROUP\n1\nRANDOM_TYPE\nuser_file\nFILE\n{grm}\n(CO)VARIANCES\n1.0\n"
            )
            if model == "h1":
                dest.write(
                    f"RANDOM_GROUP\n{total_effect+1}\nRANDOM_TYPE\nuser_file\nFILE\n{prefix}.giv\n(CO)VARIANCES\n0.5"
                )
                dest.write("\n")
            dest.write(f"OPTION method VCE")
            dest.write("\n")
            dest.write(f"OPTION maxrounds 100")
            dest.write("\n")

    def extract_logl(self, log_file):
        pattern = re.compile(r"-2logL([^0-9]+)([0-9.]+)")
        log_l = "na"
        with open(log_file) as source:
            for line in source:
                line = line.rstrip()
                if "-2logL" in line:
                    match = re.findall(pattern, line)
                    log_l = float(match[0][1])
        return log_l

    def prepare_datfile(self, phe, prefix):
        with open(prefix + ".dat", "w") as dest:
            with open(phe) as source:
                for line in source:
                    line = line.rstrip().split()
                    del line[1]
                    dest.write(" ".join(line) + "\n")

    def run_h0(self, params_file, diplo, grm, phe, model):
        prefix = grm.rstrip(".giv")
        self.prepare_params(params_file, diplo, grm, prefix, model)
        self.prepare_datfile(phe, prefix)
        command = f"mkdir {prefix} && mv {prefix}.params ./{prefix} && cp {prefix}.* {prefix} && cd {prefix} && blupf90+ {prefix}.params && cp blupf90.log ../{prefix}.log && cd .. && rm -r {prefix}"
        subprocess.call([command], shell=True)
        log_l = self.extract_logl(f"{prefix}.log")
        return log_l

    def run_blupf90(self, list_i):
        prefix, h0_logl, mid_win_point, grm = list_i
        command = f"mkdir {prefix} && mv {prefix}.params ./{prefix} && cp {prefix}.* {prefix}/ && cp {grm} {prefix}/ && cd {prefix} && blupf90+ {prefix}.params && cp blupf90.log ../{prefix}.log && cd .. && rm -r {prefix}"
        subprocess.call([command], shell=True)
        h1_logl = self.extract_logl(f"{prefix}.log")
        if h0_logl != "na" and h1_logl != "na":
            return prefix, mid_win_point, h0_logl - h1_logl


class util:

    def prepare_dat_file(self, hap_path, pheno_file, prefix, tool):
        indi_diplo_dict = {}
        with open(prefix + ".dat", "w") as dest:
            with open(hap_path) as source:
                lc = 0
                for line in source:
                    diplo = line.rstrip().split()[3]
                    indi_diplo_dict[lc] = diplo
                    lc += 1
            with open(pheno_file) as source:
                lc = 0
                for line in source:
                    line = line.rstrip().split()
                    if tool == "blupf90":
                        del line[1]
                    line.append(indi_diplo_dict[lc])
                    dest.write(" ".join(line))
                    dest.write("\n")
                    lc += 1

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
                file.write(
                    f"{i} {h1} {h2} {d} {str1} {str2}\n"
                )  # tab delimited for readability
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

    def get_maf(self, record, sample_list):
        ac = [0, 0]
        for val in sample_list:
            if record.samples[val]["GT"][0] != None:
                ac[record.samples[val]["GT"][0]] += 1
            if record.samples[val]["GT"][1] != None:
                ac[record.samples[val]["GT"][1]] += 1
        return ac


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
        outprefix,
    ):
        self.u = util()  # object containing general method
        self.asr = AsremlMethods()  # object containing asreml method
        self.blp = Blupf90Methods()  # object contaning blupf90+ method

        window_size = int(window_size)  # window size in terms of number of SNPs
        vcf = pysam.VariantFile(vcf_path)  # read vcf file
        sample_list = make_sample_list(
            pheno_file
        )  # data file of phenotypes-->sample_index, sample_id, effect1 (x1), effect2 (x2),observation (y)
        dataset = outprefix
        window_number = 0
        window_process_list = (
            []
        )  # list will collect the parameters to run asreml and blupf90
        # following condition will calculate the log-likelihood value of the null hypothesis
        if tool == "asreml":
            h0_logl = self.asr.run_h0(param_file, "na", grm, pheno_file, "h0")
        else:
            h0_logl = self.blp.run_h0(param_file, "na", grm, pheno_file, "h0")
        # Iterate through VCF records
        for i, record in enumerate(vcf):
            last_ele = -1
            is_window_size = (
                False
                if len(list(g_list_homo_dict[last_ele].keys())) < window_size
                else True
            )
            ac = self.u.get_maf(record, sample_list)
            if not record.id:
                record.id = f"{record.chrom}_{record.pos}"
            homozygosity = self.u.get_homozygosity(ac)
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
                    is_window_size = (
                        False
                        if len(list(g_list_homo_dict[last_ele].keys())) < window_size
                        else True
                    )
                else:
                    is_window_size = True
            g_list_homo_dict.append({record.id: homozygosity})
            g_list_pos_dict.append({record.id: record.pos / 1e6})
            g_list_geno_dict.append({})
            for sample in sample_list:
                sample_values = record.samples[sample]["GT"]
                g_list_geno_dict[-1][sample] = [sample_values]
            if len(g_list_homo_dict[1]) == window_size:
                sample_genotypes = g_list_geno_dict[1]
                positions = g_list_pos_dict[1]
                hzgys = g_list_homo_dict[1]
                window_number += 1
                win_min_point = np.median([list(positions.values())])
                if tool == "asreml":
                    window_process_list.append(
                        (
                            f"{dataset}.{chromosome}.{window_number}",
                            float(h0_logl),
                            win_min_point,
                        )
                    )
                else:
                    window_process_list.append(
                        (
                            f"{dataset}.{chromosome}.{window_number}",
                            float(h0_logl),
                            win_min_point,
                            grm,
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
                del g_list_geno_dict[1]
                del g_list_pos_dict[1]
                del g_list_homo_dict[1]
                if len(g_job_list) == int(num_cores):
                    with Pool(processes=len(g_job_list)) as pool:
                        pool.map(self.create_ginverse, g_job_list, 1)
                    with Pool(processes=len(window_process_list)) as pool:
                        results = pool.map(self.blp.run_blupf90, window_process_list, 1)
                    with open(outprefix + "_results.txt", "a") as dest:
                        for result in results:
                            if result:
                                dest.write(f"{chromosome} {result[1]} {result[2]}\n")
                    del g_job_list[:]
                    del window_process_list[:]
        with Pool(processes=len(g_job_list)) as pool:
            pool.map(self.create_ginverse, g_job_list, 1)
        if tool == "asreml":
            with Pool(processes=1) as pool:
                results = pool.map(self.asr.run_asreml, window_process_list, 1)
            with open(outprefix + "_results.txt", "w") as dest:
                for result in results:
                    if result:
                        dest.write(f"{chromosome} {result[1]} {result[2]}\n")
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

        self.u.prepare_dat_file(hap_path, pheno_file, prefix, tool)

        if tool == "asreml":
            self.asr.prepare_params(params_file, max_d, grm, prefix, "h1")
        else:
            self.blp.prepare_params(params_file, max_d, grm, prefix, "h1")

        command = (
            f"{dname}/cldla_snp {prefix} && {dname}/bend {prefix}.grm {prefix}.B.grm && {dname}/ginverse {max_d} {prefix}.B.grm {prefix}.giv"
            + "&& rm "
            + prefix
            + ".{Hap,Map,par,grm,B.grm}"
        )

        subprocess.call([command], shell=True)

        print(f"Generated .dat and .giv for {window_number}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="python script to create hap, map, par and giv file of windows count in cLDLA"
    )

    parser.add_argument(
        "-v", "--vcf", metavar="String", help="input phased vcf", required=True
    )
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
    parser.add_argument(
        "-p", "--pheno", metavar="String", help="phenotype file", required=True
    )
    parser.add_argument(
        "-a", "--params", metavar="String", help="parameter file", required=True
    )
    parser.add_argument(
        "-t",
        "--tool",
        default="asreml",
        const="asreml",
        nargs="?",
        choices=["asreml", "blupf90"],
        help="tools for variance component estimation --> asreml or blupf90 (default: %(default)s)",
    )

    parser.add_argument("-o", "--outprefix", help="output prefix", required=True)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
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
            args.outprefix,
        )
