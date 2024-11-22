from .reads_meta import ReadsMeta
import regex
import pandas
from collections import Counter


def hamming_dist(query, reference):
    """
    Calculate the Hamming distance between two sequences.

    :param query: The query sequence.
    :param reference: The reference sequence.
    :return: The number of positions at which the corresponding symbols are different.
    """
    mismatch_count = 0
    for _, (i, j) in enumerate(zip(query, reference)):
        if i != j:
            mismatch_count += 1
    return mismatch_count


class IdentifyAssociations(ReadsMeta):
    """
    Child class for parsing sequencing reads from LME MaPSy libraries.
    Inherits from ReadsMeta and utilizes its methods to parse reads, identify barcoded sequences,
    and match variants with mismatches.

    Inherits from:
        ReadsMeta

    Attributes
    ----------
    meta : ReadsMeta
        Instance of the ReadsMeta class containing metadata and sequencing reads.
    all_matched : pandas.DataFrame
        DataFrame of matched barcode-variant pairs.
    input_data : pandas.DataFrame
        DataFrame combining matched pairs with barcode data.

    Attributes from Parent Class (ReadsMeta)
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

    match_variants(filtered_reads, oligo_sequences, variant_seq_lengths)
        Identify barcode-variant associations with mismatches.

    """

    def __init__(self, meta):
        """
        Initialize an IdentifyAssociations instance.

        :param meta: An instance of the ReadsMeta class with metadata and sequencing reads.
        """
        _parsed_forward = meta.parse_reads(reads=meta.forward)
        _parsed_reverse = meta.parse_reads(
            reads=meta.reverse,
            reverse=True
        )

        _barcoded_sequences, self.barcodes_data = meta.parse_barcodes(
            reads=_parsed_forward,
            regex_str=meta.barcode_regex
        )

        _filtered_reverse = meta.filter_reads(
            valid_barcodes=_barcoded_sequences,
            reads=_parsed_reverse
        )

        _oligos, _variants_length = meta.parse_variants_metadata(metadata=meta.variants_info)

        # Match variants from sequencing data
        self.all_matched = self.match_variants(
            filtered_reads=_filtered_reverse,
            oligo_sequences=_oligos,
            variant_seq_lengths=_variants_length
        )

        self.input_data = self.all_matched.merge(self.barcodes_data, on='READ')

    @staticmethod
    def match_variants(
            filtered_reads: dict,
            oligo_sequences: dict,
            variant_seq_lengths: dict
    ) -> pandas.DataFrame:
        """
        Identify barcode-variant associations with mismatches.

        :param filtered_reads: Dictionary of sequencing reads containing variants.
        :param oligo_sequences: Dictionary of oligo sequences.
        :param variant_seq_lengths: Dictionary of variant lengths.
        :return: DataFrame of identified barcode-variant associations with mismatches.
        """

        matched_variants = dict()
        for read in list(filtered_reads.keys()):
            read_sequence = filtered_reads[read]

            for length in sorted(set(variant_seq_lengths.values())):
                selected_vars = {key for key, value in variant_seq_lengths.items() if value == length}
                selected_vars_sequences = {key: value for key, value in oligo_sequences.items()
                                           if key in selected_vars}
                if len(read_sequence) == 200:
                    exon_regex = regex.compile('.*(GGGATAAGACGGTAGGC){s<=5}(.{93})' +
                                               f'(.{{{length}}})' +
                                               '(.{25})(TCGTAGCACGTCACGGTTGG){s<=5}.*')
                elif len(read_sequence) != 200:
                    lenupint = len(read_sequence) - 42 - 25 - length
                    exon_regex = regex.compile('.*' + f'(.{{{lenupint - 2}}}AG)' +
                                               f'(.{{{length}}})' +
                                               '(GT.{23})(TCGTAGCACGTCACGGTTGGAGCTCCAGCCAGGTTTTCAAGC){s<=7}.*')
                match = regex.search(exon_regex, read_sequence)

                if match:
                    _, upint_match, exon_match, doint_match, _ = match.groups()
                    matched_seq = upint_match + exon_match + doint_match
                    mapped_dist = {key: hamming_dist(matched_seq, value) for key, value in
                                   selected_vars_sequences.items()}

                    mismatch_counter = Counter()
                    for mismatch in mapped_dist:
                        mismatch_counter[mapped_dist[mismatch]] += 1
                    all_generate_one = sorted({key for (key, value) in dict(mismatch_counter).items() if value == 1})

                    if len(all_generate_one) == 1:
                        first_mismatch = all_generate_one[0]
                        select_first_mismatch = \
                            list({key for (key, value) in mapped_dist.items() if value == first_mismatch})[
                                0]
                        matched_variants.update(
                            {read: [select_first_mismatch, matched_seq, first_mismatch, 0, 'SINGLE']})

                        del filtered_reads[read]
                        break

                    elif len(all_generate_one) >= 2:
                        first_mismatch, second_mismatch = all_generate_one[0], all_generate_one[1]
                        select_first_mismatch = \
                            list({key for (key, value) in mapped_dist.items() if value == first_mismatch})[
                                0]
                        if first_mismatch == 0:
                            matched_variants.update(
                                {read: [select_first_mismatch, matched_seq, first_mismatch, second_mismatch, 'SINGLE']})

                        elif first_mismatch > 0 and first_mismatch != second_mismatch:
                            matched_variants.update(
                                {read: [select_first_mismatch, matched_seq, first_mismatch, second_mismatch, 'MULTI']})

                        del filtered_reads[read]
                        break

        matched_df = pandas.DataFrame. \
            from_dict(matched_variants, orient='index'). \
            reset_index(). \
            rename(columns={
            'index': 'READ',
            0: 'VARIANT',
            1: 'MATCHED_SEQ',
            2: 'MM1',
            3: 'MM2',
            4: 'TYPE'}
        )

        return matched_df
