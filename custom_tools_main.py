import requests
import io
import os
import re
import sys
import datetime
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction as gc
from dataclasses import dataclass
from bs4 import BeautifulSoup
from dotenv import load_dotenv


def filter_fastq(input_path: str, output_filename=None,
                 gc_bounds=None, length_bounds=(0, 2 ** 32), quality_threshold=0):
    """
    Filters fastq file by GC content, length, and quality.
    :param input_path: str, path to fastq_file.
    :param output_filename: str, name of output fastq file with filtered data.
    Optional, takes input file name if not mentioned otherwise.
    :param gc_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). If no arguments are given, default value is None, returns input records.
    :param quality_threshold: float, lower limit for filtration. Default value is 0.
    :param length_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). Default value is (0, 2**32).
    :return: fastq file in fastq_filtrator_results folder.
    Raises ValueError("Too strict conditions") if no record passed the conditions.
    """

    input_path = os.path.abspath(input_path)

    if not os.path.exists('fastq_filtrator_results'):
        os.mkdir('fastq_filtrator_results')

    if output_filename is None:
        output_filename = os.path.basename(input_path)

    output_path = os.path.join('fastq_filtrator_results', output_filename)

    passed_filters = []

    for record in SeqIO.parse(input_path, "fastq"):
        seq_len = len(record.seq)
        gc_content = gc(record.seq)
        quality = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])

        # Working with GC-bounds
        if gc_bounds is not None:
            if isinstance(gc_bounds, tuple):
                if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                    continue
            else:
                if not gc_content <= gc_bounds:
                    continue

        # Working with length-bounds
        if isinstance(length_bounds, tuple):
            if not (length_bounds[0] <= seq_len <= length_bounds[1]):
                continue
        else:
            if not seq_len <= length_bounds:
                continue

        # Working with quality threshold-bounds
        if quality < quality_threshold:
            continue

        passed_filters.append(record)

    if not passed_filters:
        raise ValueError("Too strict conditions")

    with open(output_path, "w") as output_fastq:
        SeqIO.write(passed_filters, output_fastq, "fastq")


@dataclass
class GenscanOutput:
    """
    A data class for storing the output of a GENSCAN analysis.

    Attributes:
    -status (str): HTTP status code received from the GENSCAN server.
    -cds_list (dict): A dictionary containing predicted peptide sequences where keys are IDs
                                and values are sequences.
    -intron_list (dict): A dictionary of predicted introns where keys are IDs and values are the
                                start and end positions as a list of strings.
    -exon_list (dict): A dictionary of predicted exons where keys are identifiers and values are the start
                              and end positions as a list of strings.
    """
    status: str
    cds_list: dict
    intron_list: dict
    exon_list: dict

    def __repr__(self):
        cds_str = "\n".join(f"{key}: {value}" for key, value in self.cds_list.items())

        intron_str = "\n".join([f"Intron {key}: {value[0]} - {value[1]}" for key, value in self.intron_list.items()])

        exon_str = "\n".join([f"Exon {key}: {value[0]} - {value[1]}" for key, value in self.exon_list.items()])

        return (
            f"Status code: {self.status}\n\n"
            "Predicted peptides:\n"
            f"{cds_str}\n\n"
            "Predicted introns:\n"
            f"{intron_str}\n\n"
            "Predicted exons:\n"
            f"{exon_str}"
        )


def extract_exons(pre_text: str) -> dict:
    """
    Extract exon data from a <pre> HTML block returned by GENSCAN.

    Args:
    -pre_text (str): The <pre> HTML block containing the output from GENSCAN.

    Returns:
    -dict: A dictionary where keys are exon IDs and values are lists
    containing the start and end positions as strings.
    """
    exon_dict = {}

    lines = pre_text.splitlines()

    i = 0
    # find first line of the table
    while not lines[i].startswith(' 1.01'):
        i += 1

    # parse the exon table

    while not lines[i].startswith('Suboptimal exons with probability'):
        if lines[i].strip():  # skip empty lines
            data = lines[i].split()
            if data[1] in ['Init', 'Intr', 'Term', 'Sngl']:  # choose only exons
                exon_id = data[0]
                exon_begin = data[3]
                exon_end = data[4]
                exon_dict[exon_id] = [exon_begin, exon_end]
        i += 1
    return exon_dict


def calculate_introns(exons: dict) -> dict:
    """
    Calculate intron positions based on the provided exon data.

    Args:
    -exons (dict): A dictionary where keys are exon identifiers and values are lists containing
    the start and end positions as strings.

    Returns:
    -dict: A dictionary where keys are intron IDs (newly generated) and values are lists containing the
    start and end positions as strings.
    """
    introns = {}
    sorted_keys = sorted(exons, key=lambda x: (int(x.split('.')[0]), float(x)))

    last_end = None
    last_gene = None

    for key in sorted_keys:
        gene_id, exon_id = key.split('.')
        start, end = int(exons[key][0]), int(exons[key][1])

        if last_gene is None or gene_id != last_gene:
            last_gene = gene_id
            last_end = end
            continue

        if last_end is not None:
            intron_key = f"{last_gene}.{int(exon_id) - 1}"
            intron_range = [last_end + 1, start - 1]
            if intron_key not in introns:
                introns[intron_key] = [str(intron_range[0]), str(intron_range[1])]

        last_end = end

    return introns


def extract_peptides(pre_text: str) -> dict:
    """
    Extract predicted peptide sequences from a <pre> HTML text block returned by GENSCAN.

    Args:
    -pre_text (str): The <pre> HTML block containing the output from GENSCAN.

    Returns:
    -dict: A dictionary where keys are peptide IDs and values are the corresponding peptide sequences.
    """
    peptides_start = pre_text.find("Predicted peptide sequence(s):")
    peptides_end = pre_text.find("Back to GENSCAN", peptides_start)

    peptide_text = pre_text[peptides_start + len("Predicted peptide sequence(s):"):peptides_end]

    peptides = re.split(r'(?=>)', peptide_text)
    peptides = [p for p in peptides if p.strip()]

    peptide_dict = {}
    for peptide in peptides:
        header, sequence = peptide.split('\n', 1)
        header = header.split('|')[1]
        sequence = ''.join(sequence.split())
        peptide_dict[header] = sequence

    return peptide_dict


def run_genscan(sequence=None, sequence_file=None, organism="Vertebrate", exon_cutoff=1.00, sequence_name=""):
    """
    Perform a GENSCAN analysis using either an input sequence or a file containing the sequence.

    Args:
    -sequence (str, optional): The DNA sequence to analyze. Defaults to None.
    -sequence_file (str, optional): The path to a file containing the DNA sequence to analyze. Defaults to None.
    -organism (str): The type of organism the sequence is derived from. Defaults to "Vertebrate".
    Possible values are: "Vertebrate", "Arabidopsis", "Maize".
    -exon_cutoff (float): The minimum score threshold for predicted exons. Defaults to 1.00.
    Possible values are: 1.00, 0.50, 0.25, 0.10, 0.05, 0.02, 0.01.
    -sequence_name (str): An optional name for the sequence. Defaults to an empty string.

    Returns:
    -GenscanOutput: An instance of the GenscanOutput data class containing the analysis results.
    """
    base_url = "http://hollywood.mit.edu/cgi-bin/genscanw_py.cgi"
    data = {
        '-o': organism,
        '-e': exon_cutoff,
        '-n': sequence_name,
        '-p': 'Predicted peptides only'  # always
    }

    if not sequence and not sequence_file:
        raise ValueError("Either 'sequence' or 'sequence_file' must be provided.")

    if sequence:
        data['-s'] = sequence

    if sequence_file:
        files = {'-u': open(sequence_file, 'rb')}
    else:
        files = None

    response = requests.post(base_url, data=data, files=files)
    status = str(response.status_code)

    soup = BeautifulSoup(response.content, 'html.parser')

    pre_text = soup.find('pre').text

    exons = extract_exons(pre_text)
    introns = calculate_introns(exons)
    peptides = extract_peptides(pre_text)

    return GenscanOutput(status, peptides, introns, exons)


def send_telegram_message(bot_token, chat_id, message, file=None, file_name=None):
    """
        Sends a message or a document to a Telegram chat using a bot.

        Parameters:
        - chat_id (str): ID of the Telegram chat to which the message is to be sent.
        - message (str): Text message to send.
        - file (io.BytesIO, optional): File to send.
        - file_name (str, optional): Name of the file to be sent.

        Returns:
        - None
    """
    if file:
        url = f"https://api.telegram.org/bot{bot_token}/sendDocument"
        message_data = {"chat_id": chat_id, "caption": message, "parse_mode": "HTML"}
        file_data = {'document': (file_name, file)}
        requests.get(url, data=message_data, files=file_data)
    else:
        url = f"https://api.telegram.org/bot{bot_token}/sendMessage"
        message_data = {"chat_id": chat_id, "text": message, "parse_mode": "HTML"}
        requests.post(url, data=message_data)


def make_post(chat_id, bot_token, func_name, success=True, execution_time=None, output=None, exception=None):
    """
        Constructs and sends a Telegram message indicating the success or failure of a function run.

        Args:
        - chat_id (str): ID of the Telegram chat to which the message is to be sent.
        - func_name (str): Name of the function that was executed.
        - success (bool): Indicates whether the function run was successful.
        - execution_time (datetime.timedelta, optional): Time taken for the function to run.
        - output (str, optional): Output from the run of the function.
        - exception (Exception, optional): Exception object if function finished with an error.

        Returns:
        - None
        """

    post_parameters = {
        "chat_id": chat_id,
        "bot_token": bot_token,
    }

    if output:
        post_parameters["file"] = io.BytesIO(output.encode('utf-8'))
        post_parameters["file_name"] = f"{func_name}.log"

    if success:
        post_parameters["message"] = f"✅ Function <code>{func_name}</code> \
        successfully finished in <code>{execution_time}</code>."
    else:
        post_parameters["message"] = f"❌ Function <code>{func_name}</code> failed with an exception: \
        \n\n<code>{type(exception).__name__}: {str(exception)}</code>"

    send_telegram_message(**post_parameters)


def telegram_logger(chat_id):
    """
    Decorator to log function run send them to Telegram chat.
    Parameters:
    - chat_id (str): ID of the Telegram chat to which the message is to be sent.
    Returns:
    - Decorator.
    """
    load_dotenv()
    bot_token = os.getenv("TG_API_TOKEN")

    def decorator(func):
        def wrapper(*args, **kwargs):
            old_stdout = sys.stdout
            old_stderr = sys.stderr

            log_capture_string = io.StringIO()
            sys.stdout = log_capture_string
            sys.stderr = log_capture_string

            try:
                start_time = datetime.datetime.now()
                result = func(*args, **kwargs)
                end_time = datetime.datetime.now()
                output = log_capture_string.getvalue()
                execution_time = end_time - start_time

                make_post(chat_id=chat_id, bot_token=bot_token, func_name=func.__name__,
                          execution_time=execution_time, success=True, output=output)

                return result

            except Exception as e:
                output = log_capture_string.getvalue()
                make_post(chat_id=chat_id, bot_token=bot_token, func_name=func.__name__,
                          success=False, output=output, exception=e)
                raise

            finally:
                sys.stdout = old_stdout
                sys.stderr = old_stderr

        return wrapper
    return decorator
