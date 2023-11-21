# Simple script to extract the values from all the samples of a given .VCF file, assign them an id and store in a .hap file
import sys
import util
import os

# Input VCF and output hap file paths and window size (as in SNP count)
def vcf_to_custom_haplo(vcf_path, chromosome, window_size, window_number, dataset, outprefix):
    window_size = int(window_size)
    window_number = int(window_number)
    outprefix = int(outprefix)
    
    #dataset = vcf_path.split('.')[0]
    #chromosome = vcf_path.split('.')[1]

    sample_genotypes, positions, mafs, hzgys = util.read_vcf(vcf_path, window_size, window_number-1)
    
    if len(positions) != len(hzgys):
        print('Number of records mismatch')
        return exit(1) 
    
    #if window_number == 1:
    #    for i in range(int(window_size/2)):
    #        hap_path = f'{dataset}.{chromosome}.{i+1}.Hap'
    #        map_path = f'{dataset}.{chromosome}.{i+1}.Map'
    #        par_path = f'{dataset}.{chromosome}.{i+1}.par'
    # 
    #        n_samples = util.get_HAP(hap_path, sample_genotypes)
    #        util.get_MAP(map_path, positions, hzgys)
    #        util.get_PAR(par_path, window_size, window_number, n_samples)
    #
    #        print(f'Generated .hap, .map, and .par for window {i+1}')

    #else:
    #i = int(window_size/2) + window_number - 1
    hap_path = f'{dataset}.{chromosome}.{outprefix}.Hap'
    map_path = f'{dataset}.{chromosome}.{outprefix}.Map'
    par_path = f'{dataset}.{chromosome}.{outprefix}.par'

    n_samples = util.get_HAP(hap_path, sample_genotypes)
    util.get_MAP(map_path, positions, hzgys)
    util.get_PAR(par_path, window_size, outprefix, n_samples)

    print(f'Generated .hap, .map, and .par for window {window_number}')

if __name__ == "__main__":
    vcf_to_custom_haplo(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
