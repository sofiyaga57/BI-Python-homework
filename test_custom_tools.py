import unittest
import os
from sklearn.datasets import make_classification
from custom_random_forest import RandomForestClassifierCustom
from bio_files_processor import OpenFasta, FastaRecord, read_fasta_file
from custom_tools_main import run_genscan, GenscanOutput, filter_fastq

class TestRandomForestClassifierCustom(unittest.TestCase):
    """Test class for the RandomForestClassifierCustom class."""
    def test_fit_model(self):
        """Tests if the model fits with the correct number of trees."""
        X, y = make_classification(n_samples=10, random_state=42)
        model = RandomForestClassifierCustom(n_estimators=5, max_depth=3, max_features=2, random_state=42)
        model.fit(X, y)
        self.assertEqual(len(model.trees), 5)

    def test_predict_function(self):
        """Tests if predict function has the rights length of the output."""
        X, y = make_classification(n_samples=10, n_features=4, n_informative=2, n_redundant=0, random_state=42)
        X_train, y_train = X[:8], y[:8]
        X_test = X[8:]

        model = RandomForestClassifierCustom(n_estimators=5, max_depth=3, max_features=2, random_state=42)
        model.fit(X_train, y_train)
        predictions = model.predict(X_test)
        self.assertEqual(len(predictions), 2)


class TestFastaFunctionality(unittest.TestCase):
    def test_fasta_record_repr(self):
        """Test the string representation of FastaRecord."""
        record = FastaRecord(id="test_id", seq="ATCGATCGATCGAGCGCGAAAGCGCG", description="Example sequence")
        self.assertEqual(repr(record), "id='test_id', description='Example sequence...', seq='ATCGATCGATCGAGCGCGAA...'")

    def test_open_fasta_reading(self):
        """Test reading records from a FASTA file."""

        filename = "temp_test.fasta"
        with open(filename, "w") as file:
            file.write(">seq1 description1\n")
            file.write("ATGCTAGCTAGCTAGCTACA\n")
            file.write(">seq2 description2\n")
            file.write("TGCATGCTGATCGTAGCTAG\n")

        with OpenFasta(filename) as fasta_file:
            records = fasta_file.read_records()

        os.remove(filename)

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].description, "description1")
        self.assertEqual(records[0].seq, "ATGCTAGCTAGCTAGCTACA")
        self.assertEqual(records[1].id, "seq2")
        self.assertEqual(records[1].description, "description2")
        self.assertEqual(records[1].seq, "TGCATGCTGATCGTAGCTAG")


class TestRunGenscan(unittest.TestCase):
    def test_missing_sequence_and_file(self):
        """
        Test that ValueError is raised when both sequence and sequence_file parameters are None.
        """
        with self.assertRaises(ValueError) as context:
            run_genscan()
        self.assertIn("Either 'sequence' or 'sequence_file' must be provided.", str(context.exception))


class TestGenscanOutput(unittest.TestCase):
    def test_repr_method(self):
        """
        Test the __repr__ method of GenscanOutput.
        """
        output = GenscanOutput(
            status="200",
            cds_list={"ID1": "ATGCTAGCTAGCTAGCTACA", "ID2": "ATCGATCGATCGATCG"},
            intron_list={"1": ["10", "20"], "2": ["30", "40"]},
            exon_list={"A": ["1", "10"], "B": ["20", "30"]}
        )

        expected_repr = (
            "Status code: 200\n\n"
            "Predicted peptides:\n"
            "ID1: ATGCTAGCTAGCTAGCTACA\n"
            "ID2: ATCGATCGATCGATCG\n\n"
            "Predicted introns:\n"
            "Intron 1: 10 - 20\n"
            "Intron 2: 30 - 40\n\n"
            "Predicted exons:\n"
            "Exon A: 1 - 10\n"
            "Exon B: 20 - 30"
        )

        self.assertEqual(repr(output), expected_repr)


class TestReadFastaFile(unittest.TestCase):
    def test_read_fasta_file(self):
        """
        Tests if the function reads the fasta file properly.
        """
        fasta_content = """>Seq1
AGCTGCTAGCTAGCTACGATCG
>Seq2
GATCGATCGATCGTAGCTAGCTGAGCTGCTAGCTAG
"""

        with open("output_fasta.fasta", 'w') as file:
            file.write(fasta_content)

        result_dict = read_fasta_file("output_fasta.fasta")

        expected_dict = {
            'Seq1': 'AGCTGCTAGCTAGCTACGATCG',
            'Seq2': 'GATCGATCGATCGTAGCTAGCTGAGCTGCTAGCTAG'
        }

        self.assertEqual(result_dict, expected_dict)

        os.remove("output_fasta.fasta")


class TestFilterFASTQ(unittest.TestCase):
    """
    Tests if the function does proper length filtration.
    """
    def test_length_filter(self):
        with open('test_input.fastq', 'w') as f:
            f.write("@Seq1\nACGT\n+\n!!!!\n")
            f.write("@Seq2\nACGTACGT\n+\n!!!!!!!!\n")
            f.write("@Seq3\nACGTACGTACGT\n+\n!!!!!!!!!!!!\n")

        length_bounds = 8

        filter_fastq('test_input.fastq', length_bounds=length_bounds)

        with open('fastq_filtrator_results/test_input.fastq', 'r') as f:
            output_lines = f.readlines()
            num_sequences = sum(1 for line in output_lines if line.startswith('@'))

        expected_num_sequences = 2
        self.assertEqual(num_sequences, expected_num_sequences)

        os.remove('test_input.fastq')
        os.remove('fastq_filtrator_results/test_input.fastq')
        os.rmdir('fastq_filtrator_results')

if __name__ == '__main__':
    unittest.main()
