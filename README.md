# Custom bioinformatics modules


## Welcome!

Welcome to my very first bioinformatics repository! 

This repository contains the projects and homework assignments I completed during 
the one-year bioinformatics program at the Institute of Bioinformatics in Saint 
Petersburg, Russia. It includes a variety of custom bioinformatics modules 
designed to develop skills in object-oriented programming (OOP) and API usage. 

I hope you find these resources helpful and inspiring!

## Repo Structure

This repository is organized as follows to help you navigate and utilize the contents effectively:

- **Showcases.ipynb**: A Jupyter notebook that demonstrates how to use the modules listed below with practical examples.
- **custom_tools_main.py**: This script includes tools for interfacing with GenScan and Telegram through their respective APIs.
- **bio_files_processor.py**: Contains functions, classes, and a context manager designed for handling files in FASTA 
format.
- **custom_random_forest.py**: Features a custom implementation of the `RandomForestClassifier` class from scikit-learn, 
with capabilities for parallel processing.
- **test_custom_tools.py**: Contains eight simple tests to validate the functionality of the 
repository's modules.
- **data/**: Includes a folder with example data files needed for the examples demonstrated in the Jupyter notebook.

### main.py

The `main.py` module serves as the core of this repository. Below is an overview of its key functions:

- **filter_fastq**: This function filters FASTQ files based on parameters such as GC content, sequence length, and 
quality threshold.

- **run_genscan**: Performs GENSCAN analysis using either an input sequence or a sequence file.
    - class GenscanOutput: A data class designed to store the output of GENSCAN analysis, including predicted peptides, 
introns, and exons.
    - def extract_exons: Extracts exon data from HTML blocks returned by GENSCAN.
    - def calculate_introns: Calculates intron positions based on the provided exon data.
    - def extract_peptides: Extracts predicted peptide sequences from HTML text blocks returned by GENSCAN.

- **telegram_logger**: A decorator to log function runs and send them to a Telegram chat.
  - def send_telegram_message: Sends messages or documents to a Telegram chat using a bot.
  - def make_post: Constructs and sends Telegram messages indicating the success or failure of function runs.

### bio_files_processor.py

This module provides functionalities for handling FASTA and Genbank (GBK) files, including reading, writing, and 
converting them. Below is an overview of its key functions:

- **read_fasta_file**: Reads a FASTA file and saves it to a dictionary.

- **write_fasta_file**: Writes a dictionary containing FASTA data into a FASTA file.

- **convert_multiline_fasta_to_oneline**: Converts a multiline FASTA file into a one-line FASTA file.

- **select_genes_from_gbk_to_list**: Selects gene names and their translations from a Genbank (GBK) file and returns 
them as a list.

- **FastaRecord**: A data class representing a FASTA record containing sequence information.
Provides a short string representation of the FastaRecord.

- **OpenFasta**: A context manager for reading FASTA files. Supports iteration over FASTA records and reading individual
records or all records from the file.

### custom_random_forest.py

This module introduces a custom implementation of a RandomForestClassifier, capable of parallel processing for improved 
performance. 


Enjoy! ✨✨
