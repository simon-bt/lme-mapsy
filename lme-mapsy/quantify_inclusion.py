import os
import math
import pandas
import regex
import numpy
from .reads_meta import ReadsMeta


class QuantifyInclusion(ReadsMeta):
    """

    This class provides methods for quantifying the inclusion levels and filtering barcode outliers
    from LME MaPSy libraries. It extends the ReadsMeta class, which offers functionalities for
    parsing reads, identifying barcoded sequences, and matching variants with mismatches.

    Inherits from:
        ReadsMeta

    Attributes:
    -----------
    raw_inclusion_table : pandas.DataFrame
        DataFrame containing raw inclusion levels for different barcodes and variants.
    summary : pandas.DataFrame
        DataFrame summarizing inclusion counts for each barcode and variant.
    corrected_inclusion_table : pandas.DataFrame
        DataFrame containing corrected inclusion levels after filtering.
    final_summary : pandas.DataFrame
        Final summary of inclusion levels and matches after correction.

    Methods:
    --------
    __init__(meta, associations: str, condition: str, correct_aberrant: bool, replicates: bool)
        Initializes the QuantifyInclusion object with metadata and processes input data to quantify inclusion levels.

    __process_associations(cls, data_dir: str) -> dict
        Processes trusted barcode-variant associations from the given directory and returns a dictionary.

    __parse_variants_metadata(cls, metadata) -> tuple[dict, dict]
        Parses the library table for variant IDs and sequences, returning two dictionaries.

    __quantify_inclusion(cls, reads, exon_regex: str, associations: dict, trusted_barcodes: dict,
                         variants_seq: dict, variants: dict, condition: str)
                         -> tuple[pandas.DataFrame, pandas.DataFrame]
        Quantifies the inclusion levels of variants from the reads and returns the raw inclusion table and summary.

    __correct_inclusion(cls, inclusion_table, summary_table, associations)
                        -> tuple[pandas.DataFrame, pandas.DataFrame]
        Corrects inclusion levels by filtering errors and outliers, returning the corrected inclusion table and summary.
    """

    def __init__(
            self,
            meta: ReadsMeta,
            associations: str,
            condition: str,
            correct_aberrant: bool,
            replicates: bool
    ):
        """
         Initializes the QuantifyInclusion instance with metadata and processes input data to quantify inclusion levels.

        :param meta: An instance of the ReadsMeta class with metadata and sequencing reads.
        :param associations: Path to the directory containing trusted barcode-variant associations.
        :param condition: Experimental condition associated with the sequencing data.
        :param correct_aberrant: If True, perform correction for aberrant quantification.
        :param replicates: If True, consider replicates in the analysis.
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
        _variants_seq, _variants = self.__parse_variants_metadata(metadata=meta.variants_info)
        _filtered_reverse = meta.filter_reads(
            valid_barcodes=_barcoded_sequences,
            reads=_parsed_reverse
        )
        _trusted_associations = self.__process_associations(associations)
        _trusted_barcodes = {key: value for key, value in _barcoded_sequences.items() if value in _trusted_associations}
        _condition = condition
        _replicates = replicates
        self.raw_inclusion_table, self.summary = self.__quantify_inclusion(
            reads=_filtered_reverse,
            exon_regex=meta.flanking_exons_regex,
            associations=_trusted_associations,
            trusted_barcodes=_trusted_barcodes,
            variants_seq=_variants_seq,
            variants=_variants,
            condition=condition
        )
        self.corrected_inclusion_table, self.final_summary = self.__correct_inclusion(
            inclusion_table=self.raw_inclusion_table,
            summary_table=self.summary,
            associations=_trusted_associations,
            # correction=correct_aberrant
        )

    @classmethod
    def __process_associations(
            cls,
            data_dir: str
    ) -> dict:
        """
        Processes the trusted barcode-variant associations from the specified directory.

        :param data_dir: Path to the directory containing association data files.
        :return: Dictionary of trusted barcode-variant associations.
        """
        associations = dict()
        for file in os.listdir(data_dir):
            if 'TRUSTED_ASSOCIATIONS' in file:
                associations_df = pandas.read_csv(filepath_or_buffer=f"{data_dir}/{file}", sep='\t')
                for entry in associations_df.values:
                    associations.update({entry[0]: entry[1]})
        return associations

    @classmethod
    def __parse_variants_metadata(
            cls,
            metadata: str
    ) -> tuple[dict, dict]:
        """
        Parses the library table to extract variant IDs and variant sequences.

        :param metadata: Path to the library table file.
        :return: A tuple containing two dictionaries: one mapping variant sequences to IDs,
                 and another mapping variant IDs to sequences.
        """
        metadata_df = pandas.read_csv(filepath_or_buffer=metadata, sep='\t')
        variants_seqs, variants = dict(), dict()
        for entry in metadata_df[['VARIANT', 'VAR_SEQ']].values:
            variants_seqs.update({entry[1]: entry[0]})
        for entry in metadata_df[['VARIANT', 'VAR_SEQ']].values:
            variants.update({entry[0]: entry[1]})
        return variants_seqs, variants

    @classmethod
    def __quantify_inclusion(
            cls,
            reads: dict,
            exon_regex: str,
            associations: dict,
            trusted_barcodes: dict,
            variants_seq: dict,
            variants: dict,
            condition: str
    ) -> tuple[pandas.DataFrame, pandas.DataFrame]:
        """
        Quantifies the inclusion levels of variants based on the provided sequencing reads and associations.

        :param reads: Dictionary containing sequencing reads with associated variants.
        :param exon_regex: Regular expression pattern for identifying variant inclusion.
        :param associations: Dictionary of trusted barcode-variant associations.
        :param trusted_barcodes: Dictionary of trusted barcodes.
        :param variants_seq: Dictionary mapping variant sequences to IDs.
        :param variants: Dictionary mapping variant IDs to sequences.
        :param condition: Experimental condition for the analysis.
        :return: A tuple containing the raw inclusion table and an inclusion counts summary table.
        """
        trusted_reads = {key: value for key, value in reads.items() if key in trusted_barcodes}
        exc, inc = {}, {}
        match_summary = []
        for key in trusted_reads:
            variant_name = associations[trusted_barcodes[key]]
            match = regex.search(exon_regex, trusted_reads[key])

            if match:

                if len(match[0]) == 26:
                    exc[trusted_barcodes[key]] = exc.get(trusted_barcodes[key], 0) + 1

                elif len(match[0]) > 26:
                    matched_variant = match.groups()[0]

                    try:
                        _ = variants_seq[matched_variant]
                        same_seq = (matched_variant == variants[variant_name])

                        if same_seq:
                            inc[trusted_barcodes[key]] = inc.get(trusted_barcodes[key], 0) + 1
                            match_summary.append(
                                [trusted_barcodes[key], variant_name,
                                 matched_variant, variants[variant_name], same_seq])
                        else:
                            match_summary.append(
                                [trusted_barcodes[key], associations[trusted_barcodes[key]], matched_variant,
                                 'Misassigned', False])

                    except KeyError:
                        match_summary.append(
                            [trusted_barcodes[key], associations[trusted_barcodes[key]], matched_variant,
                             'Aberrant', False])

        read_counts = {**exc, **inc}
        for key, value in read_counts.items():

            if key in exc and key in inc:
                read_counts[key] = {'INC': value, 'EXC': exc[key]}

            elif key in inc and key not in exc:
                read_counts[key] = {'INC': value, 'EXC': 0}

            elif key in exc and key not in inc:
                read_counts[key] = {'INC': 0, 'EXC': value}

        counts_table = pandas.DataFrame.from_dict(read_counts, orient='index').reset_index().rename(
            columns={'index': 'BARCODE'})

        counts_table['VARIANT'] = counts_table['BARCODE'].map(associations)
        counts_table['PSI'] = (100 * counts_table['INC'] / (counts_table['INC'] + counts_table['EXC'])).round(2)

        null_variants = pandas.Series(list(set(associations.values()).difference(set(counts_table['VARIANT']))))
        null_variants_df = pandas.DataFrame(columns=['INC', 'EXC', 'PSI', 'BARCODE'], index=null_variants)
        null_variants_df.index.name = 'VARIANT'
        null_variants_df.reset_index('VARIANT', inplace=True)

        quantification_table = counts_table.append(null_variants_df)
        quantification_table['EVENT'] = quantification_table['VARIANT'].apply(lambda x: x.split('-')[1])
        quantification_table['CONDITION'] = condition
        quantification_table = quantification_table[['EVENT', 'BARCODE', 'VARIANT', 'CONDITION', 'INC', 'EXC', 'PSI']]

        summary_df = pandas.DataFrame(match_summary,
                                      columns=['BARCODE', 'VARIANT', 'MATCHED_SEQ', 'VARIANT_SEQ', 'MATCH'])
        summary_df['CONDITION'] = condition

        final_quantification_table = quantification_table.copy()
        final_summary_df = summary_df.copy()

        return final_quantification_table, final_summary_df

    @classmethod
    def __correct_inclusion(
            cls,
            inclusion_table: pandas.DataFrame,
            summary_table: pandas.DataFrame,
            associations: dict,
            # correction
    ) -> tuple[pandas.DataFrame, pandas.DataFrame]:
        """
        Corrects the inclusion levels to account for potential errors and outliers.

        :param inclusion_table: DataFrame containing per-barcode inclusion levels.
        :param summary_table: DataFrame summarizing variant inclusion matching results.
        :param associations: Dictionary of trusted barcode-variant associations.
        :return: A tuple containing the corrected inclusion table and a detailed summary table.
        """

        def __calculate_proportion(row):

            if (row['INC_ALL'] + row['ERROR']) != 0 and (row['INC_ALL'] + row['ERROR']) != numpy.nan:
                proportion = round((row['ERROR'] / (row['INC_ALL'] + row['ERROR'])) * 100, 0)

            else:
                proportion = numpy.nan
            return proportion

        def __correct_inc(row):

            if row['PROPORTION'] != numpy.nan:
                inc_corr = row['INC_TRUSTED'] * ((100 + row['PROPORTION']) / 100)

                if round(inc_corr % 1, 2) < 0.5 and inc_corr != numpy.nan:
                    inc_corr = math.floor(inc_corr)

                elif round(inc_corr % 1, 2) >= 0.5 and inc_corr != numpy.nan:
                    inc_corr = round(inc_corr)
            else:
                inc_corr = numpy.nan
            return inc_corr

        associations_df = pandas.DataFrame.from_dict(associations, orient='index'). \
            reset_index().rename(columns={'index': 'BARCODE', 0: 'VARIANT'})

        # process summary
        summaries_inc_correct = pandas.DataFrame(summary_table. \
                                                 groupby(['BARCODE', 'VARIANT', 'CONDITION'])['CONDITION'].count()). \
            rename(columns={'CONDITION': 'INC_ALL'}).reset_index()
        summaries_inc_correct['INC_ALL'] = summaries_inc_correct['INC_ALL'].astype('int64')

        # if correction:
        #     temp_errors = summary_table.query('MATCH == False'). \
        #         query('VARIANT_SEQ == \'Misassigned\'')
        #     aberrant_errors = summary_table.query('MATCH == False'). \
        #         query('VARIANT_SEQ == \'Aberrant\'')
        #     errors_final = pandas.concat([temp_errors, aberrant_errors]). \
        #         drop('VARIANT', axis=1). \
        #         merge(associations_df, on='BARCODE'). \
        #         groupby(['BARCODE', 'VARIANT', 'CONDITION'])['MATCH'].count().reset_index(). \
        #         rename(columns={'MATCH': 'ERROR'})
        #
        # else:
        errors_final = summary_table.query('MATCH == False'). \
            query('VARIANT_SEQ != \'Aberrant\''). \
            groupby(['BARCODE', 'VARIANT', 'CONDITION'])['MATCH'].count().reset_index(). \
            rename(columns={'MATCH': 'ERROR'})

        summary_inc_all = summaries_inc_correct.merge(errors_final, on=['BARCODE', 'VARIANT', 'CONDITION'],
                                                      how='outer')
        summary_inc_all['ERROR'].fillna(0, inplace=True)
        summary_inc_all['ERROR'] = summary_inc_all['ERROR'].astype('int')

        psi_table_temp = inclusion_table. \
            rename(columns={'INC': 'INC_TRUSTED'}). \
            merge(summary_inc_all, on=['BARCODE', 'VARIANT', 'CONDITION'], how='outer')

        for col_name in ['ERROR', 'INC_ALL', 'INC_TRUSTED']:
            psi_table_temp[col_name].fillna(0, inplace=True)
            psi_table_temp[col_name] = psi_table_temp[col_name].astype(int)

        psi_table_temp['PROPORTION'] = psi_table_temp.apply(__calculate_proportion, axis=1)
        psi_table_temp['PROPORTION'].fillna(0, inplace=True)
        psi_table_temp['INC_CORR'] = psi_table_temp.apply(__correct_inc, axis=1)
        psi_table_temp['INC_CORR'].fillna(0, inplace=True)
        psi_table_temp['EXC'].fillna(0, inplace=True)

        psi_table_temp['NREADS'] = psi_table_temp['INC_ALL'] + psi_table_temp['EXC']
        psi_table_temp['NREADS_CORR'] = psi_table_temp.apply(
            lambda x: x['INC_CORR'] + x['EXC'] if x['INC_ALL'] != x['ERROR'] else \
                (x['EXC'] if (x['INC_CORR'] == 0 and x['INC_ALL'] == 0 and x['ERROR'] == 0) else 0), axis=1)

        for col_name in ['NREADS', 'NREADS_CORR']:
            psi_table_temp[col_name].fillna(0, inplace=True)
            psi_table_temp[col_name] = psi_table_temp[col_name].astype('int')

        psi_table_temp['EVENT'] = psi_table_temp. \
            apply(lambda x: x['VARIANT'].split('-')[1] if pandas.notnull(x['VARIANT']) else pandas.NA, axis=1)
        return psi_table_temp, summary_inc_all
