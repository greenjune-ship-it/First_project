import unittest
import main


class MainCodeTest(unittest.TestCase):

    def test_get_nucleotides_from_fastq(self):
        allowed_symbols = ("A", "T", "G", "C", "N")
        simple_fastq = ["@S_chr13_33703600_130M/1\n",
                        "TTCCTGGTTTGTGATCCTTGGGAATATGCACCACCAAA\n",
                        "+\n",
                        "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh\n"]
        for line in main.get_nucleotides_from_fastq(simple_fastq):
            line = line.strip()
            for letter in line:
                if letter not in allowed_symbols:
                    self.fail(f"FastQ file contains unacceptable symbols")

    def test_my_grep_calc(self):
        occurence = 2
        test_pattern = "ATGCACA"
        test_lines = ("ATGCACACACA", "TGNTTTACGAA", "ATGCACAATGCACA")
        self.assertEqual(main.my_grep_calc(test_pattern, test_lines), occurence)
