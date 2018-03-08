import pandas as pd
import numpy as np
import time
import sys
import scipy
from scipy import stats
import matplotlib.pyplot as plt

class AnalyzeFeature():

    def __init__(self, similarity_matrix, featureDF):
        """Runs the electo analysis"""
        sampleRankingMat = self.getSampleRanking(similarity_matrix)
        pos_samples, neg_samples = self.getPosNeg(featureDF)

        posKS_distribution = self.KSdistribution(pos_samples,sampleRankingMat,pos_samples)
        negKS_distribution = self.KSdistribution(neg_samples,sampleRankingMat,pos_samples)

        print("posKS_distribution",posKS_distribution)
        print("negKS_distribution",negKS_distribution)

        ks_distance, ks_pvalue, direction = self.testSeparationRaw(posKS_distribution,negKS_distribution)
        print("ks_distance:", ks_distance, "ks_pvalue:", ks_pvalue, "direction:", direction)

        posHisto, negHisto = self.makeHistograms(posKS_distribution,negKS_distribution)

        # posHisto = np.array([0, 1, 1, 4, 1, 1, 0, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        # negHisto = np.array([29, 44, 28, 14,  8,  2, 11,  6,  4,  2,  0,  0,  0,  0,  0,  0,  0, 0,  0,  0])

        ############
        # plot raw #
        ############

        plt.figure(figsize=(6,4))
        plt.title("Raw KS Distribution",fontsize=16)
        plt.hist(posHisto,alpha=0.8, facecolor = 'yellow')
        plt.hist(negHisto,alpha=0.3, facecolor = 'blue')
        # plt.xlim(0,1)
        plt.savefig("/Users/Alexis/Desktop/rawKSdist.png")

        print("posHisto raw",posHisto)
        print("negHisto raw",negHisto)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        smoothNegHisto = self.smoothHistogram(negHisto)
        pseudocounts = 2
        posPseudoHisto = self.addPseudocounts(posHisto,smoothNegHisto,pseudocounts)
        smoothPosPsuedo = self.smoothHistogram(posPseudoHisto)

        print("smoothNegHisto:")
        print(smoothNegHisto)
        print("smoothPosHisto:")
        print(smoothPosPsuedo)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        #################
        # plot smoothed #
        #################
        plt.figure(figsize=(6,4))
        plt.title("Smoothed KS Distribution",fontsize=16)
        plt.hist(smoothPosPsuedo,alpha=0.8, facecolor = 'yellow')
        plt.hist(smoothNegHisto,alpha=0.3, facecolor = 'blue')
        plt.xlim(0,1)
        plt.savefig("/Users/Alexis/Desktop/smoothKSdist.png")
        #################
        """
        smooth_ks_distance, smooth_ks_pvalue, direction = self.testSeparation(pos_samples, neg_samples,pseudocounts,smoothPosPsuedo,smoothNegHisto)

        print("smooth_ks_distance",smooth_ks_distance)

        log_ratios = self.calculateLogOddsRatio(smoothNegHisto,smoothPosPsuedo)
        general_prior = len(pos_samples)/(len(neg_samples)+len(pos_samples))

        probability_vector = self.calculateProbabilities(general_prior,smoothNegHisto,smoothPosPsuedo)
        print(probability_vector)
        print("length",len(probability_vector))

        #tissue priors?
        """

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
    def getPosNeg(featureDF):
        """
        Input: Pandas DF of samples and 1's indicating positive, 0's indicating negative
        Returns list of pos and negative sample of the feature
        """
        pos_samples, neg_samples = [],[]
        D = featureDF.to_dict()

        for k,v in D.items():
            if v[0] == 0:
                neg_samples.append(k)
            elif v[0] == 1:
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

        """
        print("pos_hist",pos_hist)
        print("neg_hist",neg_hist)

        plt.figure(figsize=(6,4))
        panel1=plt.axes([0.1,0.3,.4,.333])
        panel1.hist(pos_hist,pos_bins)
        panel1.set_xlim(0,1)
        plt.title("Positive KS Distribution")

        panel2=plt.axes([0.55,0.3,.4,.666])
        panel2.hist(neg_hist, neg_bins)
        panel2.set_xlim(0,1)
        plt.title("Negative KS Distribution")
        plt.savefig("/Users/Alexis/Desktop/rawKSdist.png")
        """

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

    """
    smoothHistogram <- function(counts, breaks = c(0:20)/20, alpha = 1){

      counts_smoothed <- c()
      for(i in c(2:length(breaks))){
        #count = histogram$counts[i-1]
        bin <- breaks[i]

        bin_sum = 0
        for(j in c(2:length(breaks))){
          count <- counts[j-1]
          bin2 <- breaks[j]
          distance = abs(bin - bin2)*20

          bin_sum = bin_sum + (2^(-alpha * distance) * count)
        }
        counts_smoothed <- c(counts_smoothed,bin_sum)
      }
      counts_smoothed <- counts_smoothed/sum(counts_smoothed)
      return(counts_smoothed)
    }
    """


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

        breaks = range(21)

        for i in range(len(pos_distribution)):
            p_count = int(pos_distribution[i] * number_pos_samples)

            #getting 0 ?
            # pos_infered_data.append(runif(p_count, min=breaks[i], max=breaks[i+1]))
            p_random = np.random.uniform(low=breaks[i], high=breaks[i+1], size=p_count)
            pos_infered_data.append(p_random.tolist())

            n_count = int(round(neg_distribution[i] * number_neg_samples))
            # neg_infered_data.append(runif(n_count, min= breaks[i], max=breaks[i+1]))
            n_random = np.random.uniform(low=breaks[i], high=breaks[i+1], size=n_count)
            neg_infered_data.append(n_random.tolist())

        print("pos_infered_data")
        print(pos_infered_data)
        print("neg_infered_data")
        print(neg_infered_data)

        pos_infered_data = [item for sublist in pos_infered_data for item in sublist]
        neg_infered_data = [item for sublist in neg_infered_data for item in sublist]

        #greater of less?
        ks = scipy.stats.ks_2samp(pos_infered_data,neg_infered_data)
        ks_distance = ks[0]
        ks_pvalue = ks[1]
        direction = 1
            # ks_greater = suppressWarnings(ks.test(pos_infered_data,neg_infered_data,alternative="greater"))
            # ks_less = suppressWarnings(ks.test(pos_infered_data,neg_infered_data,alternative="less"))

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



    ############
    # ALL SNVS #
    ############

    """
    #initializing SNV df
    GBM_allSNVs = pd.DataFrame(columns = GBM_mRNA_ids)

    for line in snv_file:
        sline = line.split("\t")
        #is it the mutation we are searching for?
        # if sline[0]== mutation:
        #is it a GBM sample?
        if sline[15] in GBM_snv_ids:
            #is it an impactful type?
            if sline[8] in varianttypes:
                #add to matrix
                # GBM_ATRX[mRNAtoSNV[sline[15]]] = 1
                # s = pd.DataFrame(0, columns = GBM_mRNA_ids)
                s = pd.DataFrame(0, index=range(1), columns=GBM_mRNA_ids)
                s.name = sline[0]
                s[mRNAtoSNV[sline[15]]] = 1
                GBM_allSNVs = GBM_allSNVs.concat(s)

    print(GBM_allSNVs)
    """

    print("binary feature creation done",time.time()-startTime)

    return GBM_sim_matrix, GBM_ATRX
    # return GBM_sim_matrix, GBM_allSNVs

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
    # snv_file = pd.read_csv("/Users/Alexis/Desktop/StuartRotation/data/snvout-2.tsv",delimiter = "\t",skiprows=0, index_col=0,usecols=["Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode"])

    print("loading input files done",time.time()-startTime)


    print("Starting electo...",time.time()-startTime)

    #Running Electo on ATRX only
    GBM_sim_matrix, GBM_ATRX = GBMprep(clinical_file, IDmapfile, simMatrix, snv_file, "GBM")

    # GBM_sim_matrix, GBM_ATRX = None,None

    AnalyzeFeature(GBM_sim_matrix, GBM_ATRX)

    #Gathering all GBM SNV feature data
    # GBM_sim_matrix, GBM_allSNVs = GBMprep(clinical_file, IDmapfile, simMatrix, snv_file, "GBM")
    # AnalyzeFeature(GBM_sim_matrix, GBM_allSNVs)

    #get list of features in snv
    #for feature in feature list:
        #run AnalyzeFeature(GBM_sim_matrix, GBM_ATRX)

    print("DONE",time.time()-startTime)

if __name__ == "__main__": main()
