import pandas as pd
import numpy as np
import time
import sys
import scipy
from scipy import stats
import matplotlib.pyplot as plt
# from multiprocessing import Pool
from pathos.multiprocessing import ProcessingPool as Pool

class AnalyzeFeature():

    def runAnalyze(self, feature, similarity_matrix, featureDF):
        """Runs the electo analysis"""

        sampleRankingMat = self.getSampleRanking(similarity_matrix)
        pos_samples, neg_samples = self.getPosNeg(feature, featureDF)

        # print("feature",feature)
        # print("pos_samples",len(pos_samples))
        # print("neg_samples",len(neg_samples))


        posKS_distribution = self.KSdistribution(pos_samples,sampleRankingMat,pos_samples)
        negKS_distribution = self.KSdistribution(neg_samples,sampleRankingMat,pos_samples)

        ks_distance, ks_pvalue, direction = self.testSeparationRaw(posKS_distribution,negKS_distribution)

        posHisto, negHisto = self.makeHistograms(posKS_distribution,negKS_distribution)

        ############
        # plot raw #
        ############
        # bins = np.linspace(0,1,21) #0 to 1 with 20 bins
        # plt.figure(figsize=(6,4))
        # plt.title("Raw KS Distribution",fontsize=16)
        # plt.hist(posKS_distribution,bins,alpha=0.8, facecolor = 'yellow',label = "Positive")
        # plt.hist(negKS_distribution,bins, alpha=0.3, facecolor = 'blue',label = "Negative")
        # plt.legend(loc='upper right')
        #
        # # plt.xlim(0,1)
        # plt.savefig("/Users/Alexis/Desktop/rawKSdist.png")

        #############
        # smoothing #
        #############
        smoothNegHisto = self.smoothHistogram(negHisto)
        pseudocounts = 2
        posPseudoHisto = self.addPseudocounts(posHisto,smoothNegHisto,pseudocounts)
        smoothPosPsuedo = self.smoothHistogram(posPseudoHisto)


        #################
        # plot smoothed #
        #################
        # plt.figure(figsize=(6,4))
        # plt.title("Smoothed KS Distribution",fontsize=16)
        # plt.hist(smoothPosPsuedo,alpha=0.8, facecolor = 'yellow', label = "Positive")
        # plt.hist(smoothNegHisto,alpha=0.3, facecolor = 'blue', label = "Negative")
        # plt.legend(loc='upper right')
        # plt.xlim(0,1)
        # plt.savefig("/Users/Alexis/Desktop/smoothKSdist.png")
        #################

        smooth_ks_distance, smooth_ks_pvalue, direction = self.testSeparation(pos_samples, neg_samples,pseudocounts,smoothPosPsuedo,smoothNegHisto)
        # print("smooth",smooth_ks_distance, smooth_ks_pvalue)
        #
        #
        # print("smooth_ks_distance",smooth_ks_distance)
        # print("smooth_ks_pvalue",smooth_ks_pvalue)

        # log_ratios = self.calculateLogOddsRatio(smoothNegHisto,smoothPosPsuedo)
        # general_prior = len(pos_samples)/(len(neg_samples)+len(pos_samples))
        #
        # probability_vector = self.calculateProbabilities(general_prior,smoothNegHisto,smoothPosPsuedo)
        #
        # print("probability vector:")
        # print(probability_vector)

        return ks_distance, ks_pvalue, smooth_ks_distance, smooth_ks_pvalue

    @staticmethod
    def getSampleRanking(similarity_matrix):
        """
        calculate ranking of samples by similarity to each sample
        input: similarity matrix
        return: sample by rank matrix, containing sample names ordered by similarity to row name
        """
        samples = (list(similarity_matrix.index)) #row names
        sample_ranking_mat = pd.DataFrame(index = samples) #new matrix with same row names

        for sample in samples:

            #??? for some reason this sample breaks it ???
            if sample == 'TCGA-06-0125-01A-01R-1849-01':
                pass
            else:
                #sort row, get list of ranked sample names
                data = similarity_matrix.loc[sample]
                sorteddata = data.sort_values(ascending=False)
                L = sorteddata.index.tolist()
                sample_ranking_mat[sample] = L

            # ranked_samples = (similarity_matrix.loc[sample].sort_values(ascending=False).index.tolist())
            # sample_ranking_mat[sample] = ranked_samples

        #transform
        sample_ranking_mat = sample_ranking_mat.T

        #redo header line
        sample_ranking_mat.columns = range(len(samples))

        #remove first column
        sample_ranking_mat = sample_ranking_mat.drop(sample_ranking_mat.columns[0], axis=1)

        return sample_ranking_mat

    @staticmethod
    def getPosNeg(feature, featureDF):
        """
        Input: Pandas DF of samples and 1's indicating positive, 0's indicating negative
        Returns list of pos and negative sample of the feature
        """
        pos_samples, neg_samples = [],[]
        D = featureDF.to_dict()
        # print(D)
        for k,v in D.items():
            if v == 0:
                neg_samples.append(k)
            elif v == 1:
                pos_samples.append(k)

        return pos_samples, neg_samples

    @staticmethod
    def KStoUniform(ranked_samples, pos_samples):
        """
        calculate KS test of ranks of positive samples to uniform distribution
        input: ranked samples, names of positive samples
        return: KS distance
        """
        pos_ranks = []

        #getting enumerated dictionary of ranked samples to reference samples as rank index
        rankedsampleD = {}
        for x,y in enumerate(ranked_samples.T.values.flatten()):
            rankedsampleD[x]=y

        #getting list of positive sample rankings
        for sample in pos_samples:
            pos_ranks.append([key for key, val in rankedsampleD.items() if val == sample])

        #flattening
        pos_ranks = [val for sublist in pos_ranks for val in sublist]

        #r ks <- suppressWarnings(ks.test(pos_ranks,"punif",1,length(ranked_samples),alternative="greater"))

        # ks = scipy.stats.kstest(pos_ranks,'uniform', alternative='greater',N=len(ranked_samples))
        ks = scipy.stats.kstest(pos_ranks, 'uniform', args=(0,len(ranked_samples)-1),alternative="greater")

        return ks

    @staticmethod
    def KSdistribution(qsamples, sample_ranking_mat, pos_samples):
        """
        get KS distance for all samples in a set
        input: sample set, the sample ranking matrix, names of positive samples
        return: vector of KS distances
        """

        distances = []

        samples = (list(sample_ranking_mat.index)) #row names
        for sample in qsamples:
            if sample == 'TCGA-06-0125-01A-01R-1849-01':
                pass
            else:

                #get the rankings for this sample
                ranked_samples = sample_ranking_mat.loc[sample]
                ks = AnalyzeFeature.KStoUniform(ranked_samples, pos_samples)
                #ks[0] = ks statistic
                distances.append(ks[0])

        return distances

    @staticmethod
    def testSeparationRaw(posKS_distribution,negKS_distribution):
        """ 2 sample KS test between the positive and negative raw distributions"""

        ks = scipy.stats.ks_2samp(posKS_distribution,negKS_distribution)
        ks_distance = ks[0]
        ks_pvalue = ks[1]
        direction = 1
        return ks_distance, ks_pvalue, direction

    @staticmethod
    def makeHistograms(posKS_distribution,negKS_distribution):
        """makes histograms from distributions"""

        bins = np.linspace(0,1,21) #0 to 1 with 20 bins
        pos_hist,pos_bins=np.histogram(posKS_distribution,bins)
        neg_hist,neg_bins=np.histogram(negKS_distribution,bins)

        return pos_hist,neg_hist

    @staticmethod
    def smoothHistogram(counts):
        """smooth a histogram
        input: the histogram counts, smoothing parameter alpha (defaults to 1)
        return: smoothed histogram counts (normalized to sum up to 1)"""

        counts_smoothed = []
        alpha = 1
        breaks = np.linspace(0,1,21)

        for i in range(1,21):
            bin_ = breaks[i]
            bin_sum = 0
            for j in range(1,21):
                count = counts[j-1]
                bin2 = breaks[j]
                distance = abs(bin_ - bin2)*20

                bin_sum += (2**(-alpha*distance)*count)

            counts_smoothed.append(bin_sum)

        counts_smoothed_ = []
        for count in counts_smoothed:
            counts_smoothed_.append(count/sum(counts_smoothed))

        return counts_smoothed_

    @staticmethod
    def addPseudocounts(pos_histogram, neg_smoothed, pseudocounts):
        """ Add pseudocounts to histogram bins"""
        counts = pos_histogram.tolist()
        pos_histogram = [int(i) for i in counts]

        counts_pseudo =[]
        for i in range(len(neg_smoothed)):
            base_frequency = neg_smoothed[i]
            new_count = float(pos_histogram[i]) + (float(base_frequency) * int(pseudocounts))
            counts_pseudo.append(new_count)

        return counts_pseudo

    @staticmethod
    def testSeparation(pos_samples, neg_samples,pseudocounts,pos_distribution,neg_distribution):
        """test"""

        pos_infered_data = []
        neg_infered_data = []
        number_pos_samples = len(pos_samples) + pseudocounts
        number_neg_samples = len(neg_samples)

        breaks = np.linspace(0,1,21)

        for i in range(len(pos_distribution)):
            p_count = int(round(pos_distribution[i] * number_pos_samples))
            p_random = np.random.uniform(low=breaks[i], high=breaks[i+1], size=p_count)
            pos_infered_data.append(p_random.tolist())

            n_count = int(round(neg_distribution[i] * number_neg_samples))
            n_random = np.random.uniform(low=breaks[i], high=breaks[i+1], size=n_count)
            neg_infered_data.append(n_random.tolist())

        pos_infered_data = [item for sublist in pos_infered_data for item in sublist]
        neg_infered_data = [item for sublist in neg_infered_data for item in sublist]

        #greater of less?
        ks = scipy.stats.ks_2samp(pos_infered_data,neg_infered_data)
        ks_distance = ks[0]
        ks_pvalue = ks[1]
        direction = 1

        return ks_distance, ks_pvalue, direction

    @staticmethod
    def calculateLogOddsRatio(neg_distribution,pos_distribution):
       """calculate """
       log_ratios = []

       for i in range(1,len(neg_distribution)):
           p_neg = neg_distribution[i]
           p_pos = pos_distribution[i]
           log_ratios.append(np.log10(p_pos/p_neg))

       return log_ratios

    @staticmethod
    def calculateProbabilities(prior,neg_distribution,pos_distribution):
        """ calculate..."""
        probability_vector = []
        for i in range(len(neg_distribution)):
            p_neg = neg_distribution[i]
            p_pos = pos_distribution[i]
            prob = (p_pos*prior) / (p_pos*prior + (p_neg*(1-prior)))
            probability_vector.append(prob)

        return probability_vector

#############
# Functions #
#############
def GBMprep(clinical_file, IDmapfile, simMatrix, snv_file, disease):
    """Extract needed data from given disease and mutation"""

    ###getting ID's of GBM patients from clinical data
    GBM_patient_ids = []

    for line in clinical_file:
        sline = line.split("\t")
        #header
        if sline[0] == "bcr_patient_uuid":
            continue
        if sline[2] == disease:
            GBM_patient_ids.append(sline[1])

    ###ID mapping - getting the different IDs for the GBM samples
    GBM_mRNA_ids = []
    GBM_snv_ids = set()
    mRNAtoSNV = {}

    for line in IDmapfile:
        sline = line.split("\t")
        for sample in GBM_patient_ids:
            #check if its a GBM sample
            if sample == sline[0]:
                #must have data for both types
                if sline[6] != 'NA' and sline[3] != 'NA':
                    GBM_mRNA_ids.append(sline[6])
                    GBM_snv_ids.add(sline[3])

                    mRNAtoSNV[sline[3]] = sline[6]

    print("GBM sim matrix start",time.time()-startTime)

    ###creating GBM similarity matrix
    GBM_sim_matrix = pd.DataFrame(index=GBM_mRNA_ids, columns=GBM_mRNA_ids)

    #all needed comparisons
    comp=[(GBM_mRNA_ids[i], GBM_mRNA_ids[j]) for i in xrange(len(GBM_mRNA_ids)) for j in xrange(len(GBM_mRNA_ids))]
    for c in comp:
        GBM_sim_matrix.loc[c] = simMatrix.loc[c]

    print("GBM sim matrix done",time.time()-startTime)

    print("binary feature creation start",time.time()-startTime)

    ### Binary Feature Creation - GBM ATRX mutations
    varianttypes = ('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site')

    #############
    # ATRX only #
    ############
    """
    GBM_ATRX = pd.DataFrame(0, index = range(1), columns = GBM_mRNA_ids)
    for line in snv_file:
        sline = line.split("\t")
        #is it the mutation we are searching for?
        if sline[0]== 'ATRX':
            #is it a GBM sample?
            if sline[15] in GBM_snv_ids:
                #is it an impactful type?
                if sline[8] in varianttypes:
                    #add to matrix
                    GBM_ATRX[mRNAtoSNV[sline[15]]] = 1
    """
    # for index, row in snv_file.iterrows():
    #     if index == 'ATRX':
    #         if row['Tumor_Sample_Barcode'] in GBM_snv_ids:
    #             if row['Variant_Classification'] in varianttypes:
    #                 GBM_ATRX[mRNAtoSNV[row['Tumor_Sample_Barcode']]] = 1
    #

    ############
    # ALL SNVS #
    ############
    #initializing SNV df
    GBM_allSNVs = pd.DataFrame(columns = GBM_mRNA_ids)

    append_list = []
    for line in snv_file:
        sline = line.split("\t")
        #only GBM samples
        if sline[15] in GBM_snv_ids:
            #only proper SNV type
            if sline[8] in varianttypes:
                #add to matrix
                s = pd.DataFrame(0, index=[sline[0]], columns=GBM_mRNA_ids)
                s[mRNAtoSNV[sline[15]]] = 1
                append_list.append(s)

    GBM_allSNVs =pd.concat(append_list)

    print("collapse start",time.time()-startTime)
    #collapse the repeated rows ; set to max (which is 1)
    summ = GBM_allSNVs.groupby(GBM_allSNVs.index).max()
    print("collapse end",time.time()-startTime)

    print("filter start",time.time()-startTime)
    # get rid of genes that are not mutated or only mutation in one sample
    GBM_allSNVs = summ.loc[(summ.sum(axis=1) > 1)]
    print("filter end",time.time()-startTime)

    print("binary feature creation done",time.time()-startTime)

    # return GBM_sim_matrix, GBM_ATRX
    return GBM_sim_matrix, GBM_allSNVs

########
# MAIN #
########
def main():
    global startTime
    startTime = time.time()

    print("loading input files start",time.time()-startTime)
    clinical_file = open("/Users/Alexis/Desktop/StuartRotation/data/clinical_out.tsv","r")
    IDmapfile = open("/Users/Alexis/Desktop/StuartRotation/data/ID_mapping.tsv", "r")
    simMatrix = pd.read_csv("/Users/Alexis/Desktop/StuartRotation/data/simMatrix.pancan.atlas.imputed.tsv", delimiter = "\t",skiprows=0, index_col=0)
    snv_file = open("/Users/Alexis/Desktop/StuartRotation/data/snvout-2.tsv", "r")
    print("loading input files done",time.time()-startTime)

    #Gathering all GBM SNV feature data
    GBM_sim_matrix, GBM_allSNVs = GBMprep(clinical_file, IDmapfile, simMatrix, snv_file, "GBM")

    #get list of features in snv
    feature_list = GBM_allSNVs.index.values.tolist()

    results = pd.DataFrame(columns = ["raw_ks_distance", "raw_ks_pvalue", "smooth_ks_distance", "smooth_ks_pvalue"])

    print("Starting electo...",time.time()-startTime)

    append_list = []
    def runFeature(feature):
        feature_ = AnalyzeFeature()
        raw_ks_distance, raw_ks_pvalue, smooth_ks_distance, smooth_ks_pvalue = feature_.runAnalyze(feature, GBM_sim_matrix, GBM_allSNVs.loc[feature])
        row = pd.DataFrame(index=[feature], columns=["raw_ks_distance", "raw_ks_pvalue", "smooth_ks_distance", "smooth_ks_pvalue"])
        row["raw_ks_distance"] = raw_ks_distance
        row["raw_ks_pvalue"] = raw_ks_pvalue
        row["smooth_ks_distance"]= smooth_ks_distance
        row["smooth_ks_pvalue"] = smooth_ks_pvalue
        append_list.append(row)

    p = Pool(8)
    p.map(runFeature, feature_list)

    results = pd.concat(append_list)
    results.to_csv("/Users/Alexis/Desktop/electo_results.csv", sep=',')

    print("DONE",time.time()-startTime)

if __name__ == "__main__":
    main()
