from dataclasses import dataclass
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
import pandas
import regex


@dataclass
class ReadsMeta:
    """
     Class for handling sequencing reads and metadata.

    Attributes
    ----------
    forward : str
        Path to the forward sequencing reads file.
    reverse : str
        Path to the reverse sequencing reads file.
    sample : str
        Name of the sample.
    variants_metadata : str
        Path to the variants metadata file.
    outdir : str
        Directory for output files.
    barcode_length : int, optional
        Length of the barcode (default is 34).
    barcode_regex : str, optional
        Regex pattern for barcodes in the forward sequencing reads.
    flanking_exons_regex : str, optional
        Regex pattern for flanking exons in the reverse sequencing reads.

    Methods
    -------
    parse_reads(reads, reverse=False)
        Parse sequencing reads from a gzipped fastq file.
    parse_barcodes(reads, regex_str)
        Parse barcodes from sequencing reads using a regex pattern.
    filter_reads(valid_barcodes, reads)
        Filter reads to retain only those with valid barcodes.
    parse_variants_metadata(metadata)
        Parse variants metadata from a file.
    """

    forward: str
    reverse: str
    sample: str
    variants_metadata: str
    outdir: str
    barcode_length: int = 34
    barcode_regex: str = f".*(CTTGCTCAAC){{s<=5}}(.{{{barcode_length}}})(GAATGTCTAC){{s<=5}}.*"
    flanking_exons_regex: str = 'TGAGCGTGTTGGG(.*?)GCCAGCGAGACCG'

    @staticmethod
    def parse_reads(
            reads: str,
            reverse: bool = False
    ) -> dict:
        """
        Parse sequencing reads from a gzipped fastq file.

        :param reads: Path to the gzipped fastq file.
        :param reverse: If True, parse reverse complement of the reads (default is False).
        :return: Dictionary with read IDs as keys and sequences as values.
        """
        parsed_reads = dict()
        with gzip.open(filename=reads, mode="rt") as read:
            for title, seq, _ in FastqGeneralIterator(read):
                seq_id = title.split(' ')[0]
                if not reverse:
                    parsed_reads.update({seq_id: seq})
                else:
                    parsed_reads.update({seq_id: str(Seq(seq).reverse_complement())})
        return parsed_reads

    @staticmethod
    def parse_barcodes(
            reads: dict,
            regex_str: str
    ) -> tuple[dict, pandas.DataFrame]:
        """
       Parse barcode sequences from sequencing reads.

        :param reads: Dictionary of sequencing reads.
        :param regex_str: Regex pattern for identifying barcodes.
        :return: Tuple containing a dictionary of read IDs and their corresponding barcodes,
             and a DataFrame with columns 'READ' and 'BARCODE'.
        """
        valid_barcoded = dict()
        for key in reads:
            match = regex.search(regex_str, reads[key])
            if match:
                valid_barcoded.update({key: match.groups()[1]})
        valid_barcoded_clean = {key: value for (key, value) in valid_barcoded.items() if 'N' not in value}
        barcodes_df = pandas.DataFrame.from_dict(valid_barcoded_clean, 'index'). \
            reset_index().rename(columns={'index': 'READ', 0: 'BARCODE'})
        return valid_barcoded_clean, barcodes_df

    @staticmethod
    def filter_reads(
            valid_barcodes: dict,
            reads: dict
    ) -> dict:
        """
        Filter reads to retain only those with valid barcodes.

        :param valid_barcodes: Dictionary of valid barcoded sequences.
        :param reads: Dictionary of sequencing reads.
        :return: Dictionary of filtered sequencing reads.
        """
        keys_set = set(valid_barcodes.keys()) & set(reads.keys())
        filtered_reads = {key: reads[key] for key in keys_set}
        return filtered_reads

    @staticmethod
    def parse_variants_metadata(metadata: str) -> tuple[dict, dict]:
        """
        Parse variants metadata from a file.

        :param metadata: Path to the variants metadata file.
        :return: Tuple containing a dictionary of oligo sequences with variants as keys,
                 and a dictionary of variant lengths with variants as keys.
        """
        metadata_df = pandas.read_csv(filepath_or_buffer=metadata, sep='\t')
        oligo_sequences, variant_sequences_lengths = dict(), dict()
        for entry in metadata_df[['VARIANT', 'OLIGO']].values:
            oligo_sequences.update({entry[0]: entry[1]})
        for entry in metadata_df[['VARIANT', 'VAR_SEQ']].values:
            variant_sequences_lengths.update({entry[0]: len(entry[1])})

        return oligo_sequences, variant_sequences_lengths
