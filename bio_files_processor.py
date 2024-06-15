import os
from dataclasses import dataclass
from typing import Iterator, List


def read_fasta_file(input_fasta: str) -> dict:
    """
    Reads fasta file and saves it to dictionary
    :param input_fasta: str, fasta file
    :return dict, fasta dictionary
    """
    with open(os.path.abspath(input_fasta), mode='r') as fasta_file:
        fasta_data = {}
        current_name = None
        seq = []
        for line in fasta_file:
            line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    fasta_data[current_name] = ''.join(seq)
                current_name = line.strip('>').strip()
                seq = []
            else:
                seq.append(line.strip())
        fasta_data[current_name] = ''.join(seq)
    return fasta_data


def write_fasta_file(fasta_data: dict, output_fasta: str):
    """"
    Writes fasta dictionary into fasta file.
    :param fasta_data: dict, fasta dictionary
    :param output_fasta: str, name or path of output fasta file.
    Returns fasta file with written fasta dictionary
    """
    with open(output_fasta, mode='w') as file:
        for seq_id, seq in fasta_data.items():
            file.write('>' + seq_id + '\n')
            file.write(seq + '\n')


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta=None):
    """
    Converts multiline fasta file to oneline fasta file
    :param input_fasta: str, name of input fasta file with .fasta
    :param output_fasta: str, name of output fasta file without .fasta
    Returns file, names of sequences are seq1, seq2, etc.
    """
    input_fasta = os.path.abspath(input_fasta)
    if output_fasta is None:
        output_fasta = os.path.splitext(os.path.basename(input_fasta))[0]
    output_path = os.path.join(output_fasta + '.fasta')

    fasta_data = read_fasta_file(input_fasta)
    write_fasta_file(fasta_data=fasta_data, output_fasta=output_path)


def select_genes_from_gbk_to_list(input_gbk: str) -> list:
    """
    Selects gene's name and its translation into list
    :param input_gbk: str, name of input gbk file
    Returns list
    """

    gene_data = []
    in_cds = False
    in_translation = False
    current_gene = 'unknown'

    with open(input_gbk, mode='r') as gbk_file:
        for line in gbk_file:
            line = line.strip()
            if line.startswith("CDS"):
                in_cds = True
            if in_cds and line.startswith("/gene"):
                current_gene = line.split('gene="')[-1].split('"')[0]
            if in_cds and line.startswith("/translation"):
                in_translation = True
                current_translation = line.split('"/translation="')[-1].split('"')[1]
                continue
            if in_cds and in_translation:
                if line.endswith('"'):
                    current_translation += line.split(' ')[-1].split('"')[0]
                    in_translation = False
                    in_cds = False
                    gene_data.append(current_gene)
                    gene_data.append(current_translation)
                    current_gene = 'unknown'
                else:
                    current_translation += line

    return gene_data


@dataclass
class FastaRecord:
    """
    Represents a FASTA record containing sequence information.
    """
    id: str
    seq: str
    description: str

    def __repr__(self) -> str:
        """
        Returns short string representation of the FastaRecord.
        """
        return f"id='{self.id}', description='{self.description[:20]}...', seq='{self.seq[:20]}...'"


class OpenFasta:
    """
    Context manager for reading FASTA files.
    """
    def __init__(self, filename: str, mode: str = 'r'):
        """
        Initializes an OpenFasta instance.
        Args:
        -filename (str): The name of the FASTA file.
        -mode (str, optional): The mode in which the file is opened ('r' for reading by default).
        """
        self.filename = filename
        self.mode = mode
        self.file_handle = None
        self.current = None

    def __enter__(self) -> 'OpenFasta':
        """
        Enters the context manager and opens the FASTA file.
        Returns:
        -OpenFasta: The OpenFasta instance.
        """
        self.file_handle = open(self.filename, self.mode)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Exits the context manager and closes the file handle if open.
        Args:
        -exc_type: Exception type.
        -exc_value: Exception value.
        -traceback: Traceback information.
        """
        if self.file_handle:
            self.file_handle.close()

    def __iter__(self) -> Iterator[FastaRecord]:
        """
        Iterates over the FASTA records in the file.
        Yields:
        -FastaRecord: A FastaRecord object representing a single record.
        """
        return self

    def __next__(self) -> FastaRecord:
        """
        Retrieves the next FASTA record from the file.
        Returns:
        -FastaRecord: A FastaRecord object representing a single record.
        Raises:
        -StopIteration: If there are no more records in the file.
        """
        if self.current is None:
            line = self.file_handle.readline()
        else:
            line = self.current
            self.current = None

        if not line:
            raise StopIteration

        record_id = None
        sequence = ""
        description = ""

        while line:
            line = line.strip()
            if line.startswith(">"):
                if record_id is not None:
                    self.current = line
                    return FastaRecord(record_id, sequence, description)
                else:
                    description = ' '.join(line.split()[1:])
                    record_id = line.split()[0].strip('>')
            else:
                sequence += line

            if self.current is None:
                line = self.file_handle.readline()
            else:
                break

        if record_id:
            return FastaRecord(record_id, sequence, description)

        raise StopIteration

    def read_record(self) -> FastaRecord:
        """
        Reads the next FASTA record from the file.
        Returns:
        -FastaRecord: A FastaRecord object representing a single record.
        """
        return next(self)

    def read_records(self) -> List[FastaRecord]:
        """
        Reads all FASTA records from the file.
        Returns:
        -List[FastaRecord]: A list of FastaRecord objects representing all records in the file.
        """
        return list(self)
