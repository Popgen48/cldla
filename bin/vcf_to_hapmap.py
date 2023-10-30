# Simple script to extract the values from all the samples of a given .VCF file, assign them an id and store in a .hap file
import sys
import util
import os

# Input VCF and output hap file paths and window size (as in SNP count)
def vcf_to_custom_haplo(vcf_path, chromosome, window_size, window_number, dataset):
    window_size = int(window_size)
    window_number = int(window_number)
    
    #dataset = vcf_path.split('.')[0]
    #chromosome = vcf_path.split('.')[1]

    sample_genotypes, positions, mafs, hzgys = util.read_vcf(vcf_path, window_size, window_number-1)
    
    if len(positions) != len(hzgys):
        print('Number of records mismatch')
        return exit(1) 
    

    hap_path = f'{dataset}.{chromosome}.{window_number+1}.Hap'
    map_path = f'{dataset}.{chromosome}.{window_number+1}.Map'
    par_path = f'{dataset}.{chromosome}.{window_number+1}.par'
    
    n_samples = util.get_HAP(hap_path, sample_genotypes)
    util.get_MAP(map_path, positions, hzgys)
    util.get_PAR(par_path, window_size, window_number, n_samples)
    
    print(f'Generated .hap, .map, and .par for window {window_number+1}')

if __name__ == "__main__":
    vcf_to_custom_haplo(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
