import pandas


class ProcessIncData:
    """
    Contains methods to filter out barcode outliers and create a meta inclusion table.

    Attributes:
    -----------
    _barcode_inclusion_table : pandas.DataFrame
        DataFrame containing inclusion levels data for each barcode.
    _variants_info : str
        Path to the CSV file containing variant information.
    inclusion_table : pandas.DataFrame
        DataFrame containing filtered inclusion levels.
    accepted_barcodes : pandas.DataFrame
        DataFrame containing the barcodes accepted after filtering.
    rejected_barcodes : pandas.DataFrame
        DataFrame containing the barcodes rejected after filtering.
    meta_inclusion_table : pandas.DataFrame
        Final DataFrame combining inclusion data with variant metadata.

    Methods:
    --------
    __init__(self, barcode_inclusion_table, min_nreads, max_proportion, min_barcodes,
             min_samples, factor, variants_info)
        Initializes the ProcessIncData instance, filters barcode outliers, and prepares the meta inclusion table.

    __select_min_barcodes(cls, row, min_bc, min_cond)
        Determines if a variant should be accepted based on the minimum barcode and condition requirements.

    __filter_outliers_row(cls, row, factor)
        Determines if a barcode-variant association should be accepted based on the median PSI and PSI difference.

    __filter_outliers(cls, corr_inclusion_table, nreads, max_prop, min_barcodes,
                      min_conditions, fc)
        Filters out barcodes with aberrant inclusion levels based on provided thresholds.

    __prepare_meta_inclusion(cls, inc_table, info)
        Prepares the meta inclusion table by merging inclusion data with variant metadata.
    """

    def __init__(
            self,
            barcode_inclusion_table: pandas.DataFrame,
            min_nreads: int,
            max_proportion: float,
            min_barcodes: int,
            min_samples: int,
            factor: float,
            variants_info: str
    ):
        """
        Initialize the ProcessIncData instance.

        :param barcode_inclusion_table: DataFrame containing inclusion levels data for each barcode.
        :param min_nreads: Minimum number of reads supporting barcode-variant associations.
        :param max_proportion: Maximum proportion of corrected inclusion reads.
        :param min_barcodes: Minimum number of barcodes per variant.
        :param min_samples: Minimum number of samples containing barcode-variant association.
        :param factor: Factor for correction.
        :param variants_info: Path to the CSV file containing variant information.
        """
        self._barcode_inclusion_table = barcode_inclusion_table
        self._variants_info = variants_info
        self.inclusion_table, self.accepted_barcodes, self.rejected_barcodes = self.__filter_outliers(
            self._barcode_inclusion_table, min_nreads, max_proportion, min_barcodes, min_samples, factor
        )
        self.meta_inclusion_table = self.__prepare_meta_inclusion(self.inclusion_table, self._variants_info)

    @classmethod
    def __select_min_barcodes(
            cls,
            row: pandas.Series,
            min_bc: int,
            min_cond: int
    ) -> str:
        """
        Determines if a variant should be accepted based on the minimum barcode and condition requirements.

        :param row: DataFrame row representing a variant and its barcode counts across conditions.
        :param min_bc: Minimum number of barcodes required per variant.
        :param min_cond: Minimum number of conditions required for a barcode-variant association.
        :return: 'Accept' if the variant meets the criteria, otherwise 'Discard'.
        """
        gfp_bc = row['GFP']
        low_bc = row['LOW']
        mid_bc = row['MID']
        high_bc = row['HIGH']

        if gfp_bc >= min_bc:
            if min_cond == 3 and (low_bc >= min_bc and mid_bc >= min_bc and high_bc >= min_bc):
                call = 'Accept'
            elif min_cond == 2 and (low_bc >= min_bc and mid_bc >= min_bc) or \
                    (low_bc >= min_bc and high_bc >= min_bc) or \
                    (mid_bc >= min_bc and high_bc >= min_bc):
                call = 'Accept'
            else:
                call = 'Discard'
        else:
            call = 'Discard'
        return call

    @classmethod
    def __filter_outliers_row(
            cls,
            row: pandas.Series,
            factor: float
    ) -> str:
        """
        Determines if a barcode-variant association should be accepted based on the median PSI and PSI difference.

        :param row: DataFrame row representing a barcode-variant association with PSI values.
        :param factor: Correction factor for the PSI difference.
        :return: 'ACCEPT' if the barcode-variant association meets the criteria, otherwise 'REJECT'.
        """
        median_psi = row['MEDIAN_PSI']
        psi_diff = row['PSI_DIFF']
        if median_psi >= 90 or median_psi < 10:
            if psi_diff < 10 * factor:
                call = 'ACCEPT'
            else:
                call = 'REJECT'
        elif median_psi >= 70 or median_psi < 30:
            if psi_diff < 15 * factor:
                call = 'ACCEPT'
            else:
                call = 'REJECT'
        elif psi_diff < 20 * factor:
            call = 'ACCEPT'
        else:
            call = 'REJECT'
        return call

    @classmethod
    def __filter_outliers(
            cls,
            corr_inclusion_table: pandas.DataFrame,
            nreads: int,
            max_prop: float,
            min_barcodes: int,
            min_conditions: int,
            fc: float
    ) -> tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]:
        """
        Filters out barcodes with aberrant inclusion levels based on provided thresholds.

        :param corr_inclusion_table: DataFrame containing corrected inclusion levels.
        :param nreads: Minimum number of reads supporting barcode-variant associations.
        :param max_prop: Maximum proportion of corrected inclusion reads.
        :param min_barcodes: Minimum number of barcodes per variant.
        :param min_conditions: Minimum number of conditions containing barcode-variant association.
        :param fc: Factor for correction.
        :return: A tuple containing three DataFrames: filtered inclusion table, accepted barcodes, and rejected barcodes.
        """
        inclusion_nreads = corr_inclusion_table.query(
        f'NREADS_CORR >= {nreads} & PROPORTION <= {max_prop}')

        min_barcodes_df = inclusion_nreads[['VARIANT', 'CONDITION', 'BARCODE']]. \
            groupby(['VARIANT', 'CONDITION'])['BARCODE'].nunique().reset_index(). \
            pivot_table(index='VARIANT', columns='CONDITION', values='BARCODE').reset_index()
        min_barcodes_df['CALL'] = min_barcodes_df. \
            apply(cls.__select_min_barcodes, axis=1, args=[min_barcodes, min_conditions])
        valid_variants = min_barcodes_df[min_barcodes_df['CALL'] == 'Accept']['VARIANT'].unique()

        inclusion_filtered = inclusion_nreads[inclusion_nreads['VARIANT'].isin(valid_variants)]
        inclusion_filtered['PSI_CORR'] = round(
            inclusion_filtered['INC_CORR'] * 100 / (inclusion_filtered['INC_CORR'] + inclusion_filtered['EXC']), 2)
        inclusion_filtered['PSI_CORR'].fillna(0, inplace=True)

        psi_medians = inclusion_filtered. \
            groupby(['VARIANT', 'CONDITION'])['PSI_CORR'].median(). \
            reset_index(['VARIANT', 'CONDITION']). \
            rename(columns={'PSI_CORR': 'MEDIAN_PSI'}).round(2)
        table_merge = pandas.merge(inclusion_filtered, psi_medians, on=['VARIANT', 'CONDITION'])
        table_merge['PSI_DIFF'] = abs(table_merge['PSI_CORR'] - table_merge['MEDIAN_PSI']). \
            astype('float').round(2)
        table_merge['MEDIAN_PSI'] = table_merge['MEDIAN_PSI'].astype('float')
        table_merge['CONDITION'] = table_merge['CONDITION'].astype('category')
        table_merge['CONDITION'].cat.set_categories(['GFP', 'LOW', 'MID', 'HIGH'], inplace=True)
        table_merge.sort_values('CONDITION', inplace=True)
        table_merge['CALL'] = table_merge.apply(cls.__filter_outliers_row, axis=1, args=(fc,))

        table_merge_accepted = table_merge.query('CALL == \'ACCEPT\'')
        table_merge_rejected = table_merge.query('CALL == \'REJECT\'')

        valid_min_barcodes = table_merge_accepted[['BARCODE', 'VARIANT', 'CONDITION']]. \
            groupby(['VARIANT', 'CONDITION']).nunique().reset_index(). \
            pivot(index='VARIANT', columns='CONDITION', values='BARCODE').reset_index()
        valid_min_barcodes['CALL'] = valid_min_barcodes. \
            apply(cls.__select_min_barcodes, axis=1, args=[min_barcodes, min_conditions])
        valid_variants_filtered = valid_min_barcodes[valid_min_barcodes['CALL'] == 'Accept'][
            'VARIANT'].unique()
        valid_data = table_merge_accepted[table_merge_accepted['VARIANT'].isin(valid_variants_filtered)]

        valid_perbarcode = valid_data. \
            filter(['BARCODE', 'VARIANT', 'CONDITION', 'INC_CORR', 'EXC', 'NREADS_CORR', 'PSI_CORR'])
        valid_perbarcode['TYPE'] = valid_perbarcode['VARIANT'].apply(lambda x: 'WT' if '-WT' in x else 'VAR')
        valid_perbarcode = valid_perbarcode.set_index('BARCODE'). \
            reset_index()
        discard_index = valid_perbarcode.groupby(['VARIANT', 'CONDITION'])['BARCODE'].nunique().reset_index(). \
            query(f'BARCODE < {min_barcodes}').index
        valid_perbarcode.drop(discard_index, inplace=True)

        valid_pervar = valid_perbarcode.drop('PSI_CORR', axis=1). \
            groupby(['VARIANT', 'CONDITION']).sum(). \
            reset_index(['VARIANT', 'CONDITION'])
        valid_pervar['EVENT'] = valid_pervar['VARIANT']. \
            apply(lambda x: x.split('-')[1])
        valid_pervar['PSI'] = round(valid_pervar['INC_CORR'] * 100 / (valid_pervar['INC_CORR'] + valid_pervar['EXC']),
                                    2)
        valid_pervar['TYPE'] = valid_pervar['VARIANT'].apply(lambda x: 'WT' if '-WT' in x else 'VAR')
        # valid_pervar.drop(['INC_CORR', 'EXC', 'NREADS_CORR'], axis=1, inplace=True)  # uncomment after testing
        valid_pervar['CONDITION'] = valid_pervar['CONDITION'].astype('category')
        valid_pervar['CONDITION'].cat.set_categories(['GFP', 'LOW', 'MID', 'HIGH'], inplace=True)
        valid_pervar.sort_values(['VARIANT', 'CONDITION'], inplace=True)
        return valid_pervar, table_merge_accepted, table_merge_rejected

    @classmethod
    def __prepare_meta_inclusion(
            cls,
            inc_table: pandas.DataFrame,
            info: str
    ) -> pandas.DataFrame:
        """
        Prepares the meta inclusion table by merging inclusion data with variant metadata.

        :param inc_table: DataFrame containing per-variant inclusion levels.
        :param info: Path to the CSV file containing variant metadata.
        :return: DataFrame representing the meta inclusion table with detailed variant information.
        """
        inc_table_copy = inc_table. \
            copy()[['EVENT', 'VARIANT', 'TYPE', 'CONDITION', 'PSI', 'INC_CORR', 'EXC', 'NREADS_CORR']]
        info_df = pandas.read_csv(filepath_or_buffer=info, sep='\t'). \
            drop('OLIGO', axis=1)
        return inc_table_copy.merge(info_df, on=['VARIANT', 'EVENT'], how='left')
