import glob
import os
import re
import io
import mmap
import shutil
import argparse
import numpy as np
import pandas as pd
import multiprocessing
from functools import reduce
from time import process_time
from memory_profiler import profile


#@profile
def gunzip_fastq(file):
    import gzip
    import shutil

    with gzip.open(file, 'rb') as in_f:
        with open(file.replace('.gz', ''), 'wb') as out_f:
            shutil.copyfileobj(in_f, out_f)


#@profile
def get_nucleotides_from_fq(file_path):
    file_size = os.path.getsize(file_path)
    with open(file_path, "r+b") as f:
        map_file = mmap.mmap(f.fileno(), length=file_size, access=mmap.ACCESS_READ)
        return [line for line in iter(map_file.readline, b"")][1::4]


#@profile
def get_guides_from_file(df):
    return dict(zip(df.CODE, zip(df.GENES, df.EXONE)))


#@profile
def mmap_grep_calc(pattern, file_path):
    pattern = pattern.encode('utf-8')
    with io.open(file_path, 'r', encoding='utf-8') as f:
        file_size = os.path.getsize(file_path)
        mmap_file = mmap.mmap(f.fileno(), file_size, access=mmap.ACCESS_READ)
        return len(re.findall(pattern, mmap_file))


#@profile
def compute_results(i, nucl_f, guides_dict, tmp_folder):
    with open(f'{tmp_folder}/df_{i}.out', 'w') as out_f:
        out_f.write('gene\texone\toccurence\n')
        for guide in guides_dict:
            gene = guides_dict[guide][0]
            exone = guides_dict[guide][1]
            occurence = mmap_grep_calc(guide, nucl_f)
            result = f'{gene}\t{exone}\t{occurence}\n'
            out_f.write(result)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Search guides sequences in FASTQ file')
    parser.add_argument('-f', '--fastq', type=str, help='Path to input FAStQ file')
    parser.add_argument('-t', '--tsv', type=str, help='Path to input TSV file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    parser.add_argument('--tmp', type=str, nargs='?', default='tmp_files', help='Path to folder with temporary files')
    args = parser.parse_args()

    start_time = process_time()

    # check if file gzipped or not
    if args.fastq[-3:] == 'gz':
        gunzip_fastq(args.fastq)
        print("Input reads are unpacked.")
        args.fastq = args.fastq[:-3]
    else:
        print("Input file is recognized as unpacked.")

    # check if tmp folder exists (create or overwrite)
    if os.path.exists(args.tmp):
        print("Temporary directory exists! Overwriting...")
        shutil.rmtree(args.tmp)
        os.mkdir(args.tmp)
    else:
        os.mkdir(args.tmp)
        print("Make temporary directory.")

    print("Read FastQ file...")
    # write nucleotides in several files for reading only them in mmap
    sequences = get_nucleotides_from_fq(args.fastq)
    print("Writing nucleotides in tmp files...")
    number_of_files = 4
    sequences_split = np.array_split(sequences, number_of_files)
    for i, part in enumerate(sequences_split):
        with open(f'{args.tmp}/nucleotides_{i}.fa', 'wb') as nucl_f:
            nucl_f.write(b''.join(part))
    print("--- %s seconds for parse FastQ file ---" % (process_time() - start_time))

    guides_data = pd.read_csv(args.tsv, sep="\t")

    output_name = args.output

    dict_time = process_time()
    guides_dict = get_guides_from_file(guides_data)
    print("--- %s seconds for parse guide file ---" % (process_time() - dict_time))

    # process four separate nucleotide files and write to tmp dfs
    p = multiprocessing.Pool()
    grep_time = process_time()
    for i, nucl_f in enumerate(glob.glob(f'{args.tmp}/*.fa')):
        compute_results(i, nucl_f, guides_dict, args.tmp)
        p.apply_async(compute_results, [nucl_f])
    p.close()
    p.join()

    # summarize dfs
    with open(args.output, 'w') as out_f:
        dfs = list()
        for file in glob.glob(f'{args.tmp}/*.out'):
            # read with two index cols for suitable sum
            df = pd.read_csv(file, sep='\t', index_col=[0, 1])
            dfs.append(df)
        # string for sum multiple dfs
        dfs_sum = reduce(lambda df1, df2: df1.add(df2, fill_value=0),
                         dfs)
        dfs_sum.reset_index(level=1, inplace=True)
        dfs_sum.to_csv(out_f, sep='\t', line_terminator='\n')

    print("--- %s seconds for grep guide data ---" % (process_time() - grep_time))

    print("Remove temporary files and directory...")
    shutil.rmtree(args.tmp)

    print("--- %s seconds for whole main script ---" % (process_time() - start_time))

    print("Done. Thank you.")
