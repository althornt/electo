import pandas as pd
import time
import sys


class AnalyzeFeature():

    def __init__(self, similarity_matrix, featureDF):
        """Run electo analysis"""

        sampleRankingMat = self.getSampleRanking(similarity_matrix)

        pos_samples, neg_samples = self.getPosNeg(featureDF)

        posKS_distribution = self.KSdistribution(pos_samples)
        negKS_distribution = self.KSdistribution(neg_samples)


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


        ks = scipy.stats.kstest(pos_ranks,'norm', alternative='greater',N=len(ranked_samples))

        return ks



    def KSdistribution():
        pass

    def testSeparationRaw():
        pass

#############
# Functions #
#############
def GBMprep(clinical_file, IDmapfile, simMatrix, snv_file, disease, mutation):
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
    # IDmapfile = open("/Users/Alexis/Desktop/StuartRotation/data/ID_mapping.tsv", "r")
    GBM_ids = []
    for line in IDmapfile:
        sline = line.split("\t")
        for sample in GBM_patient_ids:
            if sample == sline[0]:
                L = [sline[0],sline[1],sline[2],sline[3],sline[4],sline[5],sline[6],sline[7],sline[8]]
                GBM_ids.append(L)

    ###getting list of GBM mRNA ids
    GBM_mRNA_ids = []
    for l in GBM_ids:
        mrna = l[6]
        if mrna == 'NA':
            pass
        else:
            GBM_mRNA_ids.append(mrna)

    ###getting list of GBM snv ids
    GBM_snv_ids = []
    for l in GBM_ids:
        snv = l[3]
        if snv == 'NA':
            pass
        else:
            GBM_snv_ids.append(snv)

    print("GBM sim matrix start",time.time()-startTime)

    #*** this part takes about 2 minutes
    ###creating GBM similarity matrix
    # mat = pd.read_csv("/Users/Alexis/Desktop/StuartRotation/data/simMatrix.pancan.atlas.imputed.tsv", delimiter = "\t",skiprows=0, index_col=False)
    # simMatrix = simMatrix.set_index('sample')
    GBM_sim_matrix = pd.DataFrame(index=GBM_mRNA_ids, columns=GBM_mRNA_ids)

    #all needed comparisons
    comp=[(GBM_mRNA_ids[i], GBM_mRNA_ids[j]) for i in xrange(len(GBM_mRNA_ids)) for j in xrange(len(GBM_mRNA_ids))]
    for c in comp:
        GBM_sim_matrix.loc[c] = simMatrix.loc[c]

    print("GBM sim matrix done",time.time()-startTime)

    print("binary feature creation start",time.time()-startTime)
    #****takes over 4 minutes
    ### Binary Feature Creation - GBM ATRX mutations
    varianttypes = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site']
    GBM_ATRX = pd.DataFrame(0, index = range(1), columns = GBM_mRNA_ids)
    # snv_file = open("/Users/Alexis/Desktop/StuartRotation/data/snvout-2.tsv", "r")

    for line in snv_file:
        sline = line.split("\t")
        if sline[0]== mutation:
            #is it a GBM sample?
            if sline[15] in GBM_snv_ids:
                #is it an impactful type?
                if sline[8] in varianttypes:
                    #add to matrix
                    GBM_ATRX[sline[15]] = 1
    print("binary feature creation done",time.time()-startTime)


    return GBM_sim_matrix, GBM_ATRX
########
# MAIN #
########
def main():
    global startTime
    startTime = time.time()

    clinical_file = open("/Users/Alexis/Desktop/StuartRotation/data/clinical_out.tsv","r")
    IDmapfile = open("/Users/Alexis/Desktop/StuartRotation/data/ID_mapping.tsv", "r")
    simMatrix = pd.read_csv("/Users/Alexis/Desktop/StuartRotation/data/simMatrix.pancan.atlas.imputed.tsv", delimiter = "\t",skiprows=0, index_col=0)
    snv_file = open("/Users/Alexis/Desktop/StuartRotation/data/snvout-2.tsv", "r")

    #Gathering GBM ATRX data
    GBM_sim_matrix, GBM_ATRX = GBMprep(clinical_file, IDmapfile, simMatrix, snv_file,"GBM", "ATRX" )

    print("Starting electo...",time.time()-startTime)

    #Running Electo
    AnalyzeFeature(GBM_sim_matrix, GBM_ATRX)

    print("DONE",time.time()-startTime)

if __name__ == "__main__": main()
