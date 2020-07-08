import argparse
import os
import sys
import numpy as np
import pandas as pd
from itertools import compress
import copy
import deTiN_utilities as du
import deTiN_SSNV_based_estimate as dssnv
import deTiN_aSCNA_based_estimate as dascna
import numpy.ma as ma

class input:
    """class which holds the required detin somatic data prior to model"""

    def __init__(self, args, ascna_probe_number_filter=200, ascna_SNP_number_filter=20, coverage_threshold=15,
                 SSNV_af_threshold=0.15, aSCNA_variance_threshold=0.025):

        # related to inputs from command line
        self.call_stats_file = args.mutation_data_path
        self.seg_file = args.cn_data_path
        self.tumor_het_file = args.tumor_het_data_path
        self.normal_het_file = args.normal_het_data_path
        self.exac_db_file = args.exac_data_path
        self.indel_file = args.indel_data_path
        self.indel_type = args.indel_data_type
        self.only_ascnas = args.only_ascnas
        if type(args.weighted_classification):
            self.weighted_classification = bool(args.weighted_classification)
        else:
            self.weighted_classification = args.weighted_classification
        if type(args.mutation_prior) == str:
            self.mutation_prior = float(args.mutation_prior)
        else:
            self.mutation_prior = args.mutation_prior
        if type(args.TiN_prior) == str:
            self.TiN_prior = float(args.TiN_prior)
        else:
            self.TiN_prior = args.TiN_prior
        if type(args.resolution) == str:
            self.resolution = int(args.resolution)
        else:
            self.resolution = args.resolution
        self.output_path = args.output_dir
        self.output_name = args.output_name
        if type(args.use_outlier_removal) == str:
            if args.use_outlier_removal.lower() == 'false':
                self.use_outlier_removal = False
            else:
                self.use_outlier_removal = True
        else:
            self.use_outlier_removal = args.use_outlier_removal
        if type(args.aSCNA_threshold) == str:
            self.aSCNA_thresh = float(args.aSCNA_threshold)
        else:
            self.aSCNA_thresh = args.aSCNA_threshold

        try:
            self.ascna_probe_number_filter = float(args.ascna_probe_number_filter)
        except AttributeError:
            self.ascna_probe_number_filter = ascna_probe_number_filter

        try:
            self.ascna_SNP_number_filter = float(args.ascna_SNP_number_filter)
        except AttributeError:
            self.ascna_SNP_number_filter = ascna_SNP_number_filter

        try:
            self.coverage_threshold = float(args.coverage_threshold)
        except AttributeError:
            self.coverage_threshold = coverage_threshold

        try:
            self.SSNV_af_threshold = float(args.SSNV_af_threshold)
        except AttributeError:
            self.SSNV_af_threshold = SSNV_af_threshold

        try:
            self.aSCNA_variance_threshold = float(args.aSCNA_variance_threshold)
        except AttributeError:
            self.aSCNA_variance_threshold = aSCNA_variance_threshold

        try:
            self.CancerHotSpotsBED = args.cancer_hot_spots
        except AttributeError:
            self.aSCNA_variance_threshold = 'NA'

        # related to inputs from class functions
        self.call_stats_table = []
        self.seg_table = []
        self.het_table = []
        self.candidates = []
        self.indel_table = []
        self.skew = 0.5

    def read_call_stats_file(self):
        fields = ['contig', 'position', 'ref_allele', 'alt_allele', 'tumor_name', 'normal_name', 't_alt_count',
                  't_ref_count'
            , 'n_alt_count', 'n_ref_count', 'failure_reasons', 'judgement','t_lod_fstar']
        fields_type = {'contig': str, 'position': np.int, 'ref_allele': str, 'alt_allele': str, 'tumor_name': str,
                       'normal_name': str,
                       't_alt_count': np.int, 't_ref_count': np.int, 'n_alt_count': np.int, 'n_ref_count': np.int,
                       'failure_reasons': str, 'judgement': str}
        try:
            self.call_stats_table = pd.read_csv(self.call_stats_file, '\t', index_col=False,
                                                comment='#', usecols=fields, dtype=fields_type)
        except (ValueError, LookupError):
            try:
                fields = ['contig', 'position', 'ref_allele', 'alt_allele', 'tumor_name', 'normal_name', 't_alt_count',
                          't_ref_count'
                    , 'n_alt_count', 'n_ref_count', 'failure_reasons', 'judgement']
                self.call_stats_table = pd.read_csv(self.call_stats_file, '\t', index_col=False,
                                                comment='#', usecols=fields, dtype=fields_type)
            except (ValueError, LookupError):
                print 'Error reading call stats skipping first two rows and trying again'
                self.call_stats_table   = pd.read_csv(self.call_stats_file, '\t', index_col=False,
                                                comment='#', skiprows=2, usecols=fields, dtype=fields_type)
        if type(self.call_stats_table['contig'][0]) == str:
            self.call_stats_table['Chromosome'] = du.chr2num(np.array(self.call_stats_table['contig']))
        else:
            self.call_stats_table['Chromosome'] = np.array(self.call_stats_table['contig']) - 1
        self.call_stats_table = self.call_stats_table[np.isfinite(self.call_stats_table['Chromosome'])]
        self.call_stats_table['genomic_coord_x'] = du.hg19_to_linear_positions(
            np.array(self.call_stats_table['Chromosome']), np.array(self.call_stats_table['position']))
        self.n_calls_in = len(self.call_stats_table)
        self.call_stats_table.reset_index(inplace=True, drop=True)

    def read_het_file(self):
        t_het_header = du.read_file_header(self.tumor_het_file)
        n_het_header = du.read_file_header(self.normal_het_file)
        cols_t_type = {t_het_header[0]: str}
        cols_n_type = {n_het_header[0]: str}
        tumor_het_table = pd.read_csv(self.tumor_het_file, '\t', index_col=False, low_memory=False, comment='#',
                                      dtype=cols_t_type)
        normal_het_table = pd.read_csv(self.normal_het_file, '\t', index_col=False, low_memory=False, comment='#',
                                       dtype=cols_n_type)
        tumor_het_table = du.fix_het_file_header(tumor_het_table)
        normal_het_table = du.fix_het_file_header(normal_het_table)
        tumor_het_table['Chromosome'] = du.chr2num(np.array(tumor_het_table['CONTIG']))
        normal_het_table['Chromosome'] = du.chr2num(np.array(normal_het_table['CONTIG']))
        tumor_het_table = tumor_het_table[np.isfinite(tumor_het_table['Chromosome'])]
        tumor_het_table['genomic_coord_x'] = du.hg19_to_linear_positions(np.array(tumor_het_table['Chromosome']),
                                                                         np.array(tumor_het_table['POSITION']))
        normal_het_table = normal_het_table[np.isfinite(normal_het_table['Chromosome'])]
        normal_het_table['genomic_coord_x'] = du.hg19_to_linear_positions(np.array(normal_het_table['Chromosome']),
                                                                          np.array(normal_het_table['POSITION']))
        tumor_het_table['AF'] = np.true_divide(tumor_het_table['ALT_COUNT'],
                                               tumor_het_table['ALT_COUNT'] + tumor_het_table['REF_COUNT'])
        normal_het_table['AF'] = np.true_divide(normal_het_table['ALT_COUNT'],
                                                normal_het_table['ALT_COUNT'] + normal_het_table['REF_COUNT'])
        self.het_table = pd.merge(normal_het_table, tumor_het_table, on='genomic_coord_x', suffixes=('_N', '_T'))

    def read_seg_file(self):
        if self.seg_file == 'NULL':
            self.seg_table = pd.DataFrame(index=[0],columns=['Chromosome','Start.bp','End.bp','n_probes','length','f','tau','genomic_coord_start','genomic_coord_end'])
            self.het_table = pd.DataFrame(index=[0],columns=['seg_id','tau','f','d','AF_T','AF_N','Chromosome','genomic_coord_x','ALT_COUNT_N'
                                                             'ALT_COUNT_T','REF_COUNT_N','REF_COUNT_T'])
        else:
            seg_header = du.read_file_header(self.seg_file)
            cols_seg_type = {seg_header[0]: str}
            self.seg_table = pd.read_csv(self.seg_file, '\t', index_col=False, low_memory=False, comment='#',
                                     dtype=cols_seg_type)
            self.seg_table = du.fix_seg_file_header(self.seg_table)
            self.seg_table['Chromosome'] = du.chr2num(np.array(self.seg_table['Chromosome']))

            self.seg_table['genomic_coord_start'] = du.hg19_to_linear_positions(np.array(self.seg_table['Chromosome']),
                                                                            np.array(self.seg_table['Start.bp']))
            self.seg_table['genomic_coord_end'] = du.hg19_to_linear_positions(np.array(self.seg_table['Chromosome']),
                                                                          np.array(self.seg_table['End.bp']))

    def annotate_call_stats_with_allelic_cn_data(self):
        f_acs = np.zeros([self.n_calls_in, 1]) + 0.5
        tau = np.zeros([self.n_calls_in, 1]) + 2
        for i, r in self.seg_table.iterrows():
            f_acs[np.logical_and(np.array(self.call_stats_table['genomic_coord_x']) >= r['genomic_coord_start'],
                                 np.array(self.call_stats_table['genomic_coord_x']) <= r['genomic_coord_end'])] = r.f
            tau[np.logical_and(np.array(self.call_stats_table['genomic_coord_x']) >= r['genomic_coord_start'],
                               np.array(self.call_stats_table['genomic_coord_x']) <= r[
                                   'genomic_coord_end'])] = r.tau + 0.001
        self.call_stats_table['tau'] = tau
        self.call_stats_table['f_acs'] = f_acs

    def annotate_het_table(self):
        seg_id = np.zeros([len(self.het_table), 1]) - 1
        tau = np.zeros([len(self.het_table), 1]) + 2
        f = np.zeros([len(self.het_table), 1]) + 0.5
        for seg_index, seg in self.seg_table.iterrows():
            het_index = np.logical_and(self.het_table['genomic_coord_x'] >= seg['genomic_coord_start'],
                                       self.het_table['genomic_coord_x'] <= seg['genomic_coord_end'])
            ix = list(compress(xrange(len(het_index)), het_index))
            seg_id[ix] = seg_index
            tau[ix] = seg['tau']
            f[ix] = seg['f']
        self.het_table['seg_id'] = seg_id
        self.het_table['tau'] = tau
        self.het_table['f'] = f
        d = np.ones([len(self.het_table), 1])
        d[np.array(self.het_table['AF_T'] <= 0.5, dtype=bool)] = -1
        self.skew = 0.5
        self.het_table['d'] = d
    def read_and_preprocess_SSNVs(self):
        self.read_call_stats_file()
        self.read_seg_file()
        self.annotate_call_stats_with_allelic_cn_data()
        if not self.indel_file == 'None':
            if not self.indel_type == 'None':
                self.indel_table = du.read_indel_vcf(self.indel_file, self.seg_table, self.indel_type)
            else:
                print 'Warning: if indels are provided you must also specify indel data source using --indel_data_type'
                print 'no indels will be returned'
                self.indel_file = 'None'
                self.indel_type = 'None'
    def read_and_preprocess_aSCNAs(self):
        self.read_seg_file()
        self.read_het_file()
        self.seg_table = du.filter_segments_based_on_size_f_and_tau(self.seg_table, self.aSCNA_thresh,
                                                                    self.ascna_probe_number_filter)
        self.annotate_het_table()
        self.het_table = du.remove_sites_near_centromere_and_telomeres(self.het_table)

    def read_and_preprocess_data(self):
        self.read_and_preprocess_SSNVs()
        self.read_and_preprocess_aSCNAs()




class output:
    """ combined from deTiN's models
    reclassified SSNVs based on TiN estimate are labeled KEEP in judgement column
    self.SSNVs['judgement'] ==  KEEP

    confidence intervals (CI_tin_high/low) represent 95% interval
    """

    def __init__(self, input, ssnv_based_model, ascna_based_model):

        # previous results
        self.input = input
        self.ssnv_based_model = ssnv_based_model
        self.ascna_based_model = ascna_based_model
        # useful outputs
        self.SSNVs = input.candidates
        self.joint_log_likelihood = np.zeros([self.input.resolution, 1])
        self.joint_posterior = np.zeros([self.input.resolution, 1])
        self.CI_tin_high = []
        self.CI_tin_low = []
        self.TiN = []
        self.p_null = 1
        # variables
        self.TiN_range = np.linspace(0, 1, num=self.input.resolution)
        self.TiN_int = 0
        # threshold for accepting variants based on the predicted somatic assignment
        # if p(S|TiN) exceeds threshold we keep the variant.
        self.threshold = 0.5
        # defines whether to remove events based on predicted exceeding predicted allele fractions
        # if Beta_cdf(predicted_normal_af;n_alt_count+1,n_ref_count+1) <= 0.01 we remove the variant
        self.use_outlier_threshold = input.use_outlier_removal
        if self.input.indel_table.isnull().values.sum() == 0:
            self.indels = self.input.indel_table

    def calculate_joint_estimate(self):
        # do not use SSNV based estimate if it exceeds 0.3 (this estimate can be unreliable at high TiNs due to
        # germline events)
        if self.ssnv_based_model.TiN <= 0.3 and ~np.isnan(self.ascna_based_model.TiN):
            if len(self.ascna_based_model.centroids) > 1:
                reselect_cluster = np.argmin(np.abs(self.ascna_based_model.centroids / 100 - self.ssnv_based_model.TiN))
                self.ascna_based_model.TiN_likelihood = self.ascna_based_model.cluster_TiN_likelihoods[reselect_cluster]
                print('reselected cluster based on SSNVs')
            # combine independent likelihoods
            self.joint_log_likelihood = self.ascna_based_model.TiN_likelihood + self.ssnv_based_model.TiN_likelihood
            # normalize likelihood to calculate posterior
            self.joint_posterior = np.exp(self.ascna_based_model.TiN_likelihood + self.ssnv_based_model.TiN_likelihood
                                          - np.nanmax(
                self.ascna_based_model.TiN_likelihood + self.ssnv_based_model.TiN_likelihood))
            self.joint_posterior = np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))
            self.CI_tin_low = self.TiN_range[next(x[0] for x in enumerate(
                np.cumsum(np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                                                  x[1] > 0.025)]
            self.CI_tin_high = self.TiN_range[
                next(x[0] for x in enumerate(np.cumsum(
                    np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                     x[1] > 0.975)]

            self.TiN_int = np.nanargmax(self.joint_posterior)
            self.TiN = self.TiN_range[self.TiN_int]

            zero_tin_ssnv_model = copy.deepcopy(self.ssnv_based_model)
            zero_tin_ssnv_model.TiN = 0
            zero_tin_ssnv_model.expectation_of_z_given_TiN()
            zero_tin_ssnv_model.maximize_TiN_likelihood()
            zero_total_l = zero_tin_ssnv_model.TiN_likelihood + self.ascna_based_model.TiN_likelihood
            zero_total_l = np.exp(zero_total_l - np.nanmax(zero_total_l))
            self.p_null = np.true_divide(zero_total_l,np.nansum(zero_total_l))[0]
            print 'joint TiN estimate = ' + str(self.TiN)
        # use only ssnv based model
        elif ~np.isnan(self.ascna_based_model.TiN):
            # otherwise TiN estimate is = to aSCNA estimate
            print 'SSNV based TiN estimate exceed 0.3 using only aSCNA based estimate'
            self.joint_log_likelihood = self.ascna_based_model.TiN_likelihood
            self.joint_posterior = np.exp(
                self.ascna_based_model.TiN_likelihood - np.nanmax(self.ascna_based_model.TiN_likelihood))
            self.joint_posterior = np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))
            self.CI_tin_low = self.TiN_range[next(x[0] for x in enumerate(
                np.cumsum(np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                                                  x[1] > 0.025)]
            self.CI_tin_high = self.TiN_range[
                next(x[0] for x in enumerate(np.cumsum(
                    np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                     x[1] > 0.975)]
            self.TiN_int = np.nanargmax(self.joint_posterior)
            self.TiN = self.TiN_range[self.TiN_int]
            self.p_null = self.joint_posterior[0]
        # use only aSCNA based estimate
        elif ~np.isnan(self.ssnv_based_model.TiN) and self.ssnv_based_model.TiN <= 0.3:
            print 'No aSCNAs only using SSNV based model'
            self.joint_log_likelihood = self.ssnv_based_model.TiN_likelihood
            self.joint_posterior = np.exp(
                self.ssnv_based_model.TiN_likelihood - np.nanmax(self.ssnv_based_model.TiN_likelihood))
            self.joint_posterior = np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))
            self.CI_tin_low = self.TiN_range[next(x[0] for x in enumerate(
                np.cumsum(np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                                                  x[1] > 0.025)]
            self.CI_tin_high = self.TiN_range[
                next(x[0] for x in enumerate(np.cumsum(
                    np.ma.masked_array(np.true_divide(self.joint_posterior, np.nansum(self.joint_posterior))))) if
                     x[1] > 0.975)]
            self.TiN_int = np.nanargmax(self.joint_posterior)
            self.TiN = self.TiN_range[self.TiN_int]
            zero_tin_ssnv_model = copy.deepcopy(self.ssnv_based_model)
            zero_tin_ssnv_model.TiN = 0
            zero_tin_ssnv_model.expectation_of_z_given_TiN()
            zero_tin_ssnv_model.maximize_TiN_likelihood()
            zero_total_l = zero_tin_ssnv_model.TiN_likelihood
            zero_total_l = np.exp(zero_total_l - np.nanmax(zero_total_l))
            self.p_null = np.true_divide(zero_total_l, np.nansum(zero_total_l))[0]
        else:
            print 'insuffcient data to generate TiN estimate.'
            self.CI_tin_high = 0
            self.CI_tin_low = 0
            self.joint_posterior = np.zeros([self.input.resolution, 1])
            self.joint_posterior[0] = 1
            self.TiN_int = 0
            self.TiN = 0
            self.p_null = 1
        pH1 = self.joint_posterior[self.TiN_int]
        #print(self.joint_posterior)
        #print(self.p_null)
        # code to deal with underflows
        if ma.is_masked(self.p_null):
            self.p_null = 0
        pH0 = self.p_null

        p_model = np.true_divide(self.input.TiN_prior * pH1,
                          (self.input.TiN_prior * pH1) + ((1 - self.input.TiN_prior) * pH0))
        if p_model < 0.5 or ~np.isfinite(p_model):
            print 'insufficient evidence to justify TiN > 0'
            self.joint_posterior = np.zeros([self.input.resolution, 1])
            self.joint_posterior[0] = 1
            self.TiN_int = 0
            self.TiN = 0
            self.CI_tin_high = 0
            self.CI_tin_low = 0

    def reclassify_mutations(self):
        # calculate p(Somatic | given joint TiN estimate)
        if self.input.weighted_classification == True:
            numerator = np.zeros(len(self.ssnv_based_model.p_TiN_given_S))
            denominator = np.zeros(len(self.ssnv_based_model.p_TiN_given_S))
            for idx, p in enumerate(self.joint_posterior):
                if p > 0.001:
                    num_iter = (p * self.ssnv_based_model.p_somatic * self.ssnv_based_model.p_TiN_given_S[:, idx])
                    numerator = numerator + num_iter
                    denom_iter = num_iter + p * (np.array(
                        [1 - self.ssnv_based_model.p_somatic] * np.nan_to_num(
                            self.ssnv_based_model.p_TiN_given_G[:, idx])))
                    denominator = denominator + denom_iter

        else:
            numerator = self.ssnv_based_model.p_somatic * np.expand_dims(
                self.ssnv_based_model.p_TiN_given_S[:, self.TiN_int], 1)
            denominator = numerator + np.array(
                [1 - self.ssnv_based_model.p_somatic] * np.expand_dims(
                    np.nan_to_num(self.ssnv_based_model.p_TiN_given_G[:, self.TiN_int]), 1))
        self.SSNVs.loc[:, ('p_somatic_given_TiN')] = np.nan_to_num(np.true_divide(numerator, denominator))
        # expected normal allele fraction given TiN and tau
        af_n_given_TiN = np.multiply(self.ssnv_based_model.tumor_f, self.ssnv_based_model.CN_ratio[:, self.TiN_int])
        # probability of normal allele fraction less than or equal to predicted fraction
        self.SSNVs.loc[:, 'p_outlier'] = self.ssnv_based_model.rv_normal_af.cdf(af_n_given_TiN + 0.01)
        if self.TiN_int == 0:
            print 'Estimated 0 TiN no SSNVs will be recovered outputing deTiN statistics for each site'
        elif self.use_outlier_threshold:
            # remove outliers mutations p(af_n >= E[af_n|TiN]) < 0.05
            self.SSNVs['judgement'][np.logical_and(self.SSNVs['p_somatic_given_TiN'] > self.threshold,
                                                   self.SSNVs['p_outlier'] >= 0.01)] = 'KEEP'
        else:
            self.SSNVs['judgement'][self.SSNVs['p_somatic_given_TiN'] > self.threshold] = 'KEEP'
        if self.input.indel_file != 'None':
            if self.input.indel_table.isnull().values.sum() == 0:
                indel_model = dssnv.model(self.input.indel_table, self.input.mutation_prior, self.input.resolution)
                indel_model.generate_conditional_ps()
                self.indels = self.input.indel_table
                numerator = indel_model.p_somatic * np.expand_dims(indel_model.p_TiN_given_S[:, self.TiN_int], 1)
                denominator = numerator + np.array(
                    [1 - indel_model.p_somatic] * np.expand_dims(np.nan_to_num(
                        indel_model.p_TiN_given_G[:, self.TiN_int]), 1))
                af_n_given_TiN = np.multiply(indel_model.tumor_f, indel_model.CN_ratio[:, self.TiN_int])
                self.indels.loc[:, ('p_somatic_given_TiN')] = np.nan_to_num(np.true_divide(numerator, denominator))
                self.indels.loc[:, 'p_outlier'] = indel_model.rv_normal_af.cdf(af_n_given_TiN)
                if self.TiN_int == 0:
                    print 'Estimated 0 TiN no indels will be recovered outputing deTiN statistics for each site'
                elif self.use_outlier_threshold:
                    # remove outliers mutations p(af_n >= E[af_n|TiN]) < 0.05
                    self.indels['filter'][np.logical_and(self.indels['p_somatic_given_TiN'] > self.threshold,
                                                     self.indels['p_outlier'] >= 0.01)] = 'PASS'
                else:
                    self.indels['filter'][self.indels['p_somatic_given_TiN'] > self.threshold] = 'PASS'
            elif self.input.indel_table.isnull().values.sum() >  0:
                self.indels = self.input.indel_table

__version__ = '1.0'


def main():
    """ deTiN pipeline. Method operates in two stages (1) estimating tumor in normal via candidate SSNVs and SCNAS.
        (2) Performing variant re-classification using bayes rule.
    """

    parser = argparse.ArgumentParser(description='Estimate tumor in normal (TiN) using putative somatic'
                                                 ' events see Taylor-Weiner & Stewart et al. 2017')
    # input files
    parser.add_argument('--mutation_data_path',
                        help='Path to mutation candidate SSNV data.'
                             'Supported formats: MuTect call-stats', required=False, default = 'NULL')
    parser.add_argument('--cn_data_path',
                        help='Path to copy number data.'
                             'Supported format: AllelicCapseg .seg file. Generated by GATK4 AllelicCNV.',
                        required=False, default='NULL')
    parser.add_argument('--tumor_het_data_path',
                        help='Path to heterozygous site allele count data in tumor. Generated by GATK4 GetBayesianHetCoverage.'
                             'Required columns: CONTIG,POS,REF_COUNT and ALT_COUNT', required=False,
                        default='NULL')
    parser.add_argument('--normal_het_data_path',
                        help='Path to heterozygous site allele count data in normal. Generated by GATK4 GetBayesianHetCoverage'
                             'Required columns: CONTIG,POS,REF_COUNT and ALT_COUNT', required=False,
                        default = 'NULL')
    parser.add_argument('--exac_data_path',
                        help='Path to exac af > 0.01 pickle. Can be generated by downloading ExAC VCF and running build_exac_pickle',
                        required=False)
    parser.add_argument('--indel_data_path',
                        help='Path to candidate indels data.'
                             'Supported formats: Strelka / MuTect2 VCFs', required=False, default='None')
    parser.add_argument('--indel_data_type',
                        help='MuTect2 or Strelka'
                             'Caller used to generate indels', required=False, default='None')
    # output related arguments
    parser.add_argument('--output_name', required=True,
                        help='sample name')
    parser.add_argument('--output_dir', help='directory to put plots and TiN solution', required=False, default='.')
    # model related parameters
    parser.add_argument('--mutation_prior', help='prior expected ratio of somatic mutations to rare germline events'
                        , required=False, default=0.05)
    parser.add_argument('--aSCNA_threshold', help='minor allele fraction threshold for calling aSCNAs.'
                        , required=False, default=0.1)
    parser.add_argument('--TiN_prior',
                        help='expected frequency of TiN contamination in sequencing setting used for model selection',
                        required=False, default=0.5)
    parser.add_argument('--use_outlier_removal',
                        help='remove sites from recovered SSNVs where allele fractions significantly exceed predicted fraction',
                        required=False, default=True)
    parser.add_argument('--resolution',
                        help='number of TiN bins to consider default = 101 corresponds to 0.01 TiN levels'
                        , required=False, default=101)
    parser.add_argument('--weighted_classification',
                        help='integrate variant classification over all values of TiN'
                        , required=False, default=False)
    parser.add_argument('--ascna_probe_number_filter', help='number of probes to require for an aSCNA to be considered'
                        , required=False, default=200)
    parser.add_argument('--ascna_SNP_number_filter', help='number of probes to require for an aSCNA to be considered'
                        , required=False, default=20)
    parser.add_argument('--coverage_threshold', help='number of reads required to use a site for TiN estimation',
                        required=False,
                        default=15)
    parser.add_argument('--SSNV_af_threshold', help='fraction of alternate alleles required for site to be used '
                                                    'for SSNV TiN estimation',
                        required=False,
                        default=.2)
    parser.add_argument('--aSCNA_variance_threshold',
                        help='variance of segment allele shift tolerated before removing segment '
                             'as artifact', required=False, default=0.025)
    parser.add_argument('--cancer_hot_spots',
                        help='Optional BED file of cancer hot spot mutations which the user has a stronger prior on being somatic e.g. BRAF v600E mutations.'
                             'The format of this file is Chromosome\tPosition\tProbability. Note this will override the mutation prior at these locations'
                        , required=False, default='NA')
    parser.add_argument('--only_ascnas',
                        help='only use ascna data for TiN estimation',required=False, action='store_true')
    args = parser.parse_args()
    if args.cn_data_path == 'NULL' and args.mutation_data_path == 'NULL':
        print 'One of CN data or SSNV data are required.'
        sys.exit()
    elif args.cn_data_path =='NULL':
        di = input(args)
        di.read_and_preprocess_SSNVs()
        di.candidates = du.select_candidate_mutations(di.call_stats_table, di.exac_db_file)
        ssnv_based_model = dssnv.model(di.candidates, di.mutation_prior, di.resolution, di.SSNV_af_threshold,
                                       di.coverage_threshold, di.CancerHotSpotsBED, skew = di.skew)
        ssnv_based_model.perform_inference()
        ascna_based_model = dascna.model(di.seg_table, di.het_table, di.resolution)
        ascna_based_model.TiN = np.nan

    elif args.mutation_data_path=='NULL':
        di = input(args)
        di.read_and_preprocess_aSCNAs()
        di.candidates = pd.DataFrame(index=[0],columns=['contig', 'position', 'ref_allele', 'alt_allele', 'tumor_name', 'normal_name',
                't_alt_count','t_ref_count', 'n_alt_count', 'n_ref_count', 'failure_reasons', 'judgement','genomic_coord_x','f_acs','tau'])
        ssnv_based_model = dssnv.model(di.candidates, di.mutation_prior, di.resolution, di.SSNV_af_threshold,
                                       di.coverage_threshold, di.CancerHotSpotsBED)
        ssnv_based_model.TiN = np.nan
        ascna = False
        # identify aSCNAs and filter hets
        if len(di.seg_table) > 0:
            di.aSCNA_hets = du.ensure_balanced_hets(di.seg_table, di.het_table)
            if len(di.aSCNA_hets) > 0:
                di.aSCNA_segs, di.convergent_segs = du.identify_aSCNAs(di.seg_table, di.aSCNA_hets, di.aSCNA_thresh,
                                                                       di.ascna_SNP_number_filter,
                                                                       di.aSCNA_variance_threshold)
                if len(di.aSCNA_segs) > 0:
                    ascna_based_model = dascna.model(di.aSCNA_segs, di.aSCNA_hets, di.resolution)
                    ascna_based_model.perform_inference()
                    ascna = True
        if not ascna:
            ascna_based_model = dascna.model(di.seg_table, di.het_table, di.resolution)
            ascna_based_model.TiN = np.nan
    else:
        di = input(args)
        di.read_and_preprocess_data()
        # identify candidate mutations based on MuTect flags.
        # kept sites are flagged as KEEP or rejected for normal lod and/or alt_allele_in_normal
        di.candidates = du.select_candidate_mutations(di.call_stats_table, di.exac_db_file)
        n_calls_pre = np.sum(di.candidates['judgement'] == "KEEP")
        # generate SSNV based model using candidate sites
        ssnv_based_model = dssnv.model(di.candidates, di.mutation_prior, di.resolution, di.SSNV_af_threshold,
                                   di.coverage_threshold, di.CancerHotSpotsBED)
        ssnv_based_model.perform_inference()
        if di.only_ascnas == True:
            ssnv_based_model.TiN = np.nan
            print 'Only using aSCNA data'
        ascna = False
        # identify aSCNAs and filter hets
        if len(di.seg_table) > 0:
            di.aSCNA_hets = du.ensure_balanced_hets(di.seg_table, di.het_table)
            if len(di.aSCNA_hets) > 0:
                di.aSCNA_segs,di.convergent_segs = du.identify_aSCNAs(di.seg_table, di.aSCNA_hets, di.aSCNA_thresh, di.ascna_SNP_number_filter,
                                               di.aSCNA_variance_threshold)
                if len(di.aSCNA_segs) > 0:
                    ascna_based_model = dascna.model(di.aSCNA_segs, di.aSCNA_hets, di.resolution)
                    ascna_based_model.perform_inference()
                    ascna = True
        if not ascna:
            ascna_based_model = dascna.model(di.seg_table, di.het_table, di.resolution)
            ascna_based_model.TiN = np.nan

        # combine models and reclassify mutations
    do = output(di, ssnv_based_model, ascna_based_model)
    do.calculate_joint_estimate()
    if len(do.SSNVs)>1:
        do.reclassify_mutations()
        do.SSNVs.drop('Chromosome', axis=1, inplace=True)
    n_calls_post = np.sum(do.SSNVs['judgement']=="KEEP")
    n_calls_added = n_calls_post - n_calls_pre
    # make output directory if needed
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # write deTiN reclassified SSNVs
    do.SSNVs.to_csv(path_or_buf=do.input.output_path + '/' + do.input.output_name + '.deTiN_SSNVs.txt', sep='\t',
                    index=None)

    if not di.indel_file == 'None':
        #if 'Chromosome' in do.indels.columns:
        do.indels.drop('Chromosome', axis=1, inplace=True)
        do.indels.to_csv(path_or_buf=do.input.output_path + '/' + do.input.output_name + '.deTiN_indels.txt', sep='\t',
                         index=None)
    # write plots
    if not np.isnan(ascna_based_model.TiN):
        do.ascna_based_model.segs['Chromosome'] = do.ascna_based_model.segs['Chromosome'] + 1
        do.ascna_based_model.segs.to_csv(path_or_buf=do.input.output_path + '/' + do.input.output_name + '.deTiN_aSCNAs.txt', sep='\t',
                    index=None)
        du.plot_kmeans_info(ascna_based_model, do.input.output_path, do.input.output_name)
        du.plot_TiN_models(do)
        du.plot_aSCNA_het_data(do)
    if not np.isnan(ssnv_based_model.TiN):
        du.plot_SSNVs(do)
    # write TiN and CIs
    file = open(do.input.output_path + '/' + do.input.output_name + '.TiN_estimate.txt', 'w')
    file.write('%s' % (do.TiN))
    file.close()

    file = open(do.input.output_path + '/' + do.input.output_name + '.TiN_estimate_CI.txt', 'w')
    file.write('%s - %s' % (str(do.CI_tin_low), str(do.CI_tin_high)))
    file.close()

    file = open(do.input.output_path + '/' + do.input.output_name + '.number_of_SSNVs_added.txt','w')
    file.write('%s\n'% int(n_calls_added))

if __name__ == "__main__":
    main()
