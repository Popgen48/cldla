import multiprocessing
import sys
from vcf_to_hapmap import vcf_to_custom_haplo

with open('record_counts.txt', 'r') as text_file:
    line = text_file.readline()
    record_count = int(line.split(':')[1].strip()) if line.split(':')[0].strip() == sys.argv[1] else None

if record_count is None:
    print('Record count not found. Please check \'record_counts.txt\'')
    exit(1)

def parallelize(vcf_path, window_size, n_proc):

    arg_list = [(vcf_path, window_size, window_number) for window_number in range(record_count - int(window_size) + 1)]
    print(f'Total number of windows to process: {len(arg_list)}')

    pool = multiprocessing.Pool(int(n_proc))
    pool.map(vcf_to_custom_haplo, arg_list)
    pool.close()
    pool.join()
    
    print('All processes finished')
    
if __name__ == '__main__':
    parallelize(sys.argv[1], sys.argv[2], sys.argv[3])