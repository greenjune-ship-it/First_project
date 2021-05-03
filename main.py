import argparse
import re
import pandas as pd


def get_nucleotides_from_fastq(lines):
    nucleotide_lines = lines[1::4]
    return nucleotide_lines


def get_guides_from_file(df):
    return dict(zip(df.CODE, zip(df.GENES, df.EXONE)))


def my_grep_calc(pattern, lines):
    return len(tuple(filter(lambda line: re.search(pattern, line), lines)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Search guides sequences in FASTQ file')
    parser.add_argument('-f', '--fastq', type=str, help='Path to input FAStQ file')
    parser.add_argument('-t', '--tsv', type=str, help='Path to input TSV file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    args = parser.parse_args()

    with open(args.fastq) as in_f:
        fastq_lines = in_f.readlines()
        sequences = get_nucleotides_from_fastq(fastq_lines)

    guides_data = pd.read_csv(args.tsv, sep="\t")
    output_name = args.output

    guides_dict = get_guides_from_file(guides_data)

    with open(output_name, "w") as out_f:
        for guide in guides_dict:
            gene = guides_dict[guide][0]
            exone = guides_dict[guide][1]
            result_line = f"{gene}\t{exone}\t{my_grep_calc(guide, sequences)}\n"
            out_f.write(result_line)
    print("Done.")
