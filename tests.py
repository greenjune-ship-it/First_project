import unittest
import main


class MainCodeTest(unittest.TestCase):

    def test_get_nucleotides_from_fastq(self):
        allowed_symbols = ("A", "T", "G", "C", "N")
        simple_fastq = 'test_data/test_get_nucleotides_from_fastq.txt'
        bytes_str = main.get_nucleotides_from_fq(simple_fastq)
        for line in bytes_str:
            line = line.decode("utf-8").strip()
            for letter in line:
                if letter not in allowed_symbols:
                    self.fail("FastQ file contains unacceptable symbols")

    def test_my_grep_calc_1(self):
        occurence =
        simple_nucl_file = 'test_data/test_my_grep_calc_1.txt'
        test_pattern = "ATGCACA"
        self.assertEqual(main.mmap_grep_calc(test_pattern, simple_nucl_file), occurence)
