from pathlib import Path
import pandas
from .reads_meta import ReadsMeta

def calc_match_pct(
        row: pandas.Series,
        max_mismatch_pct: int
) -> str:
    """
    Calculate the match percentage for a given row of data.

    :param row: DataFrame row containing mismatch data.
    :param max_mismatch_pct: Maximum mismatch percentage allowed.
    :return: 'ACCEPT' or 'REJECT' based on the specified criteria.
    """

    mm1 = row['MM1']
    mm2 = row['MM2']
    match_type = row['TYPE']
    match_len = len(row['MATCHED_SEQ'])
    pct_mismatched = (mm1 * 100) / match_len

    if match_type == 'SINGLE' and (mm1 == 0 or mm2 == 0) and pct_mismatched <= max_mismatch_pct:
        call = 'ACCEPT'
    elif match_type == 'MULTI' and mm1 == 1 and mm2 == 3 and pct_mismatched <= max_mismatch_pct:
        call = 'ACCEPT'
    elif mm2 > 2 and mm2 / mm1 >= 2 and pct_mismatched <= max_mismatch_pct:
        call = 'ACCEPT'
    else:
        call = 'REJECT'
    return call

class ResolveAssociations:
    """
    Class to handle the resolution of barcode-variant associations.

    This class processes barcode-variant association data, filters based on mismatch percentages,
    and resolves associations while classifying them as accepted or rejected.

    Attributes:
    -----------
    accepted_mismatches : pd.DataFrame
        DataFrame storing accepted mismatched associations.
    rejected_mismatches : pd.DataFrame
        DataFrame storing rejected mismatched associations.
    correct_associations : pd.DataFrame
        DataFrame storing correct barcode-variant associations.
    misassignments : pd.DataFrame
        DataFrame storing misassigned barcode-variant associations.
    input_trusted_associations : pd.DataFrame
        DataFrame storing the path to the trusted associations input file.
    resolved_associations : pd.DataFrame
        DataFrame storing resolved barcode-variant associations.

    Methods:
    --------
    __init__(input_dir: str, max_mismatch_pct: int, min_nreads: int, nreads_second_highest: int,
             pct_second_best: int, min_rep: int)
        Initializes the ResolveAssociations object with input directory and parameters for filtering and resolving associations.

    call_mismatches(data: pd.DataFrame, pct_mms: int) -> tuple[pd.DataFrame, pd.DataFrame]
        Static method to call mismatches on the input data.

    __read_input_data(cls, associations_dir: str) -> pd.DataFrame
        Reads the input trusted associations data from files in the specified directory.

    __filter_associations(cls, accepted_data: pd.DataFrame, min_nreads: int) -> tuple[pd.DataFrame, pd.DataFrame]
        Filters associations based on the minimum number of reads and returns correct assignments and misassignments.

    __resolve_associations(cls, correct_data: pd.DataFrame, misassignment_data: pd.DataFrame,
                           nreads_second_highest: int, pct_second_best: int, min_rep: int) -> tuple[pd.DataFrame, pd.DataFrame]
        Resolves barcode-variant associations based on several criteria and returns resolved and trusted associations.
    """

    def __init__(
            self,
            input_dir: str,
            max_mismatch_pct: int,
            min_nreads: int,
            nreads_second_highest: int,
            pct_second_best: int,
            min_rep: int
    ):
        """
        Initialize the ResolveAssociations instance.
        This method initializes the object with the specified parameters and processes the input data to
        determine the accepted and rejected mismatched associations, correct assignments, misassignments,
        and resolves barcode-variant associations.

        :param input_dir: Directory containing the input trusted associations data files.
        :param max_mismatch_pct: Maximum mismatch percentage allowed for accepting an association.
        :param min_nreads: Minimum number of reads required for an association to be considered.
        :param nreads_second_highest: Number of reads of the second most supported variants allowed.
        :param pct_second_best: Maximum percentage of reads contributing to the second most supported variants.
        :param min_rep: Minimum number of replicates with a barcode-variant association.
        """

        _input_data: pandas.DataFrame = self.__read_input_data(associations_dir=input_dir)
        _max_mismatch_pct: int = max_mismatch_pct
        _min_nreads: int = min_nreads
        _nreads_second_highest: int = nreads_second_highest
        _pct_second_best: int = pct_second_best
        _min_rep: int = min_rep

        self.accepted_mismatches, self.rejected_mismatches = self.call_mismatches(
            data=_input_data,
            pct_mms=_max_mismatch_pct
        )
        self.correct_associations, self.misassignments = self.__filter_associations(
            accepted_data=self.accepted_mismatches,
            min_nreads=_min_nreads
        )
        self.input_trusted_associations, self.resolved_associations = self.__resolve_associations(
            correct_data=self.correct_associations,
            misassignment_data=self.misassignments,
            nreads_second_highest=_nreads_second_highest,
            pct_second_best=_pct_second_best,
            min_rep=_min_rep
        )


    @staticmethod
    def call_mismatches(
            data: pandas.DataFrame,
            pct_mms: int
    ) -> tuple[pandas.DataFrame, pandas.DataFrame]:
        """
        Determine accepted and rejected mismatches in the input data.

        :param data: Concatenated input DataFrame containing barcode-variant associations.
        :param pct_mms: Maximum allowed percentage of mismatches.
        :return: A tuple containing two DataFrames: one with accepted associations and another with rejected associations.
        """
        data_copy = data.copy()
        data_copy['CALL'] = data_copy.apply(
            calc_match_pct,
            axis=1,
            args=[pct_mms]
        )
        accepted_associations = data_copy.query('CALL == \'ACCEPT\'')
        rejected_associations = data_copy.query('CALL == \'REJECT\'')
        return accepted_associations, rejected_associations

    @classmethod
    def __read_input_data(
            cls,
            associations_dir: str
    ) -> pandas.DataFrame:
        """
        Reads all input data required for processing barcode-variant associations.

        :param associations_dir: Path to the directory containing association data files.
        :return: Concatenated DataFrame of barcode-variant associations with sample information.
        """

        files = [f for f in Path(associations_dir).iterdir() if f.suffix == '.tab']
        dfs = []
        for file in files:
            sample_name = file.stem.split(sep='_', maxsplit=3)[-1]
            df_temp = pandas.read_csv(file, sep='\t')
            df_temp['SAMPLE'] = sample_name
            dfs.append(df_temp)
        return pandas.concat(dfs)

    @classmethod
    def __filter_associations(
            cls,
            accepted_data: pandas.DataFrame,
            min_nreads: int
    ) -> tuple[pandas.DataFrame, pandas.DataFrame]:
        """
        Filters barcode-variant associations based on read counts and uniqueness.

        :param accepted_data: DataFrame of barcode-variant associations.
        :param min_nreads: Minimum number of reads required to retain an association.
        :return: DataFrame of filtered barcode-variant associations.
        """
        accepted_data_copy = accepted_data.copy()
        filtered_nreads = accepted_data_copy.groupby(['BARCODE', 'VARIANT', 'SAMPLE'])['READ'].count(). \
            reset_index().query(f'READ >= {min_nreads}')

        correct_assignments = filtered_nreads.groupby(['BARCODE', 'SAMPLE'])['VARIANT'].nunique(). \
            reset_index(). \
            query('VARIANT == 1'). \
            drop('VARIANT', axis=1). \
            merge(filtered_nreads, on=['BARCODE', 'SAMPLE'], how='left')

        misassignments = filtered_nreads.groupby(['BARCODE', 'SAMPLE'])['VARIANT'].nunique(). \
            reset_index().query('VARIANT > 1').drop('VARIANT', axis=1). \
            merge(filtered_nreads, on=['BARCODE', 'SAMPLE'], how='left')
        return correct_assignments, misassignments

    @classmethod
    def __resolve_associations(
            cls,
            correct_data: pandas.DataFrame,
            misassignment_data: pandas.DataFrame,
            nreads_second_highest: int,
            pct_second_best: int,
            min_rep: int
    ) -> tuple[pandas.DataFrame, pandas.DataFrame]:
        """
        Resolves barcode-variant misassignments and merges dataframes as needed to ensure correct associations.

        :param correct_data: DataFrame of correct barcode-variant associations without misassignments.
        :param misassignment_data: DataFrame of barcode-variant associations with misassignments.
        :param nreads_second_highest: Number of reads of the second most supported variants.
        :param pct_second_best: Maximum percentage of reads contributing to the second most supported variants.
        :param min_rep: Minimum number of replicates with a barcode-variant association.
        :return: A tuple containing two DataFrames: trusted associations and resolved associations.
        """
        correct_data_copy = correct_data.copy()
        misassignment_data_copy = misassignment_data.copy()

        temp_filtered_misassignments = []
        for sample in misassignment_data_copy['SAMPLE'].unique():
            df_sample = misassignment_data_copy.query(f'SAMPLE == \'{sample}\'')
            keep_barcodes = []

            for barcode in df_sample['BARCODE'].unique():
                barcode_df = df_sample[df_sample['BARCODE'] == barcode]
                second_largest_support = list(barcode_df['READ'].nlargest(2))[1]

                if (barcode_df.shape[0] > 2) & (second_largest_support <= nreads_second_highest):
                    largest_df = barcode_df[barcode_df['READ'] == barcode_df['READ'].nlargest(1).values[0]]
                    keep_barcodes.append(largest_df)

                elif barcode_df.shape[0] == 2:
                    max_nreads = barcode_df['READ'].max()
                    percent = round((second_largest_support * 100) / max_nreads + second_largest_support, 2)

                    if percent <= pct_second_best:
                        largest_df = barcode_df[barcode_df['READ'] == barcode_df['READ'].nlargest(1).values[0]]
                        keep_barcodes.append(largest_df)

            if not len(keep_barcodes) == 0:
                keep_df = pandas.concat(keep_barcodes)
                temp_filtered_misassignments.append(keep_df)

        if not len(temp_filtered_misassignments) == 0:
            filtered_misassignments = pandas.concat(temp_filtered_misassignments)
            all_assignments = pandas.concat([correct_data_copy, filtered_misassignments])
        else:
            all_assignments = correct_data_copy.copy()

        assignments_min_rep = all_assignments. \
            groupby(['BARCODE'])['SAMPLE']. \
            nunique(). \
            reset_index(). \
            query(f'SAMPLE >= {min_rep}'). \
            drop('SAMPLE', axis=1). \
            merge(all_assignments, on=['BARCODE'], how='left')

        resolved_assignments = assignments_min_rep. \
            groupby(['BARCODE'])['VARIANT']. \
            nunique(). \
            reset_index().query('VARIANT == 1'). \
            drop('VARIANT', axis=1). \
            merge(assignments_min_rep, on='BARCODE', how='left')

        trusted_associations = resolved_assignments[['BARCODE', 'VARIANT']]
        trusted_associations = trusted_associations[~trusted_associations.duplicated()]
        return trusted_associations, resolved_assignments
