import sys
import os.path
import re
from multiprocessing import Pool
import subprocess
from vcf_to_local_lrt import AsremlMethods, Blupf90Methods


class MatToLrt:

    def __init__(self):
        self.asmthod = AsremlMethods()
        self.blpmthod = Blupf90Methods()

    def main_func(self, cpu, csv, param, pheno, prefix, tool):
        grm_line = True
        chrom, grm = "na", "na"
        h1_inputs = []
        h0_logl = "na"
        self.cpu = int(cpu)
        with open(csv) as source:
            for line in source:
                line = line.rstrip().split(",")
                if grm_line:
                    chrom, grm = line
                    grm_line = False
                    h0_logl = self.run_h0(chrom, grm, param, pheno, tool)
                else:
                    mid_win_point, win_prefix, drm, dat, max_diplo = line
                    new_prefix = self.prepare_h1_params(
                        chrom,
                        dat,
                        drm,
                        grm,
                        max_diplo,
                        param,
                        pheno,
                        prefix,
                        tool,
                    )
                    h1_inputs.append([new_prefix, h0_logl, mid_win_point, grm])
        if tool == "asreml":
            results = self.run_h1_asreml(h1_inputs)
        if tool == "blupf90":
            results = self.run_h1_blupf90(h1_inputs)
        with open(f"{prefix}.{chrom}.results.txt", "a") as dest:
            for result in results:
                if result:
                    dest.write(f"{chrom} {result[0]} {result[1]} {result[2]}\n")

    def run_h1_asreml(self, list_i):
        with Pool(processes=1) as pool:
            results = pool.map(self.asmthod.run_asreml, list_i, 1)
        return results

    def run_h1_blupf90(self, list_i):
        with Pool(processes=self.cpu) as pool:
            results = pool.map(self.blpmthod.run_blupf90, list_i, 1)
        return results

    def run_h0(self, chrom, grm, param, pheno, tool):
        command = f"ln -s {grm}"
        subprocess.call([command], shell=True)
        grm = os.path.basename(grm)
        log_l = "na"
        if tool == "asreml":
            log_l = self.asmthod.run_h0(param, 999, grm, pheno, "h0")
        if tool == "blupf90":
            log_l = self.blpmthod.run_h0(param, 999, grm, pheno, "h0")
        if log_l == "na":
            print("logl for H0 cannot be 0")
            sys.exit(1)
        return log_l

    def prepare_h1_params(self, chrom, dat, drm, grm, max_diplo, param, pheno, prefix, tool):
        drm_fn = os.path.basename(drm)
        pattern = re.compile(r"(.*)\.([0-9]+)\.giv")
        match = re.findall(pattern, drm_fn)
        new_prefix = f"{prefix}.{chrom}.{match[0][1]}"
        command = f"cp {drm} {new_prefix}.giv"
        subprocess.call([command], shell=True)
        self.create_new_pheno(pheno, new_prefix, dat)
        dat = os.path.basename(dat)
        grm = os.path.basename(grm)
        if tool == "blupf90":
            self.blpmthod.prepare_params(param, max_diplo, grm, new_prefix, "h1")
        if tool == "asreml":
            self.asmthod.prepare_params(param, max_diplo, grm, new_prefix, "h1")
        return new_prefix

    def create_new_pheno(self, new_pheno, new_prefix, old_pheno):
        iid_diplo = {}
        header = True
        with open(old_pheno) as source:
            for line in source:
                line = line.rstrip().split()
                iid_diplo[line[0]] = line[-1]
        with open(f"{new_prefix}.dat", "w") as dest:
            header = True
            with open(new_pheno) as source:
                for line in source:
                    if header:
                        header = False
                        pass
                    else:
                        line = line.rstrip().split()
                        line.append(iid_diplo[line[0]])
                        dest.write("\t".join(line) + "\n")


if __name__ == "__main__":
    mtl = MatToLrt()
    mtl.main_func(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
