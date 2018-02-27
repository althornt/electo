import pandas as pd
import time


class AnalyzeFeature():

    def __init__(self,featureName, similarity_matrix, featureDF):
        """Run electo analysis"""

        featureName.sampleRankingMat = self.getSampleRanking(featureDF)

        featureName.pos_samples, featureName.neg_samples = self.get_pos_neg(featureDF)

        featureName.posKS_distribution = self.KSdistribution(featureName.pos_samples)
        featureName.negKS_distribution = self.KSdistribution(featureName.neg_samples)


    @staticmethod
    def getSampleRanking(similarity_matrix):
        """
        calculate ranking of samples by similarity to each sample
        input: similarity matrix
        return: sample by rank matrix, containing sample names ordered by similarity to row name
        """
        samples = (list(similarity_matrix.index)) #row names
        sample_ranking_mat = pd.DataFrame(index = samples) #new matrix with same row names

        print("samples length", len(samples))

        for sample in samples:
            #sort row, get list of ranked sample names
            ranked_samples = (similarity_matrix.loc[sample].sort_values(ascending=False).index.tolist())
            sample_ranking_mat[sample] = ranked_samples

        #transform
        sample_ranking_mat = sample_ranking_mat.T

        #redo header line
        sample_ranking_mat.columns = range(len(samples))

        #remove first column
        sample_ranking_mat = sample_ranking_mat.drop(sample_ranking_mat.columns[0], axis=1)

        return sample_ranking_mat


    def getPosNeg(featureDF):
        """
        Input: Pandas DF of samples and 1's indicating positive, 0's indicating negative
        Returns list of pos and negative sample of the feature
        """

        pos_samples, neg_samples = [],[]
        D = GBM_ATRX.to_dict()

        for k,v in D.items():
            if v[0] == 0:
                neg_samples.append(k)
            elif v[0] == 1:
                pos_samples.append(k)

        return pos_samples, neg_samples


    def KStoUniform():
        pass

    def KSdistribution():
        pass

    def testSeparationRaw():
        pass

def GBMprep(clinical, IDmapfile, mat, snv_file):
    ##################
    # Data wrangling #
    ##################

    ###getting ID's of GBM patients from clinical data
    # clinical = open("/Users/Alexis/Desktop/StuartRotation/data/clinical_out.tsv","r")
    GBM_patient_ids = []

    for line in clinical:
        sline = line.split("\t")
        #header
        if sline[0] == "bcr_patient_uuid":
            continue
        if sline[2] == "GBM":
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

    #*** this part takes about 2 minutes
    ###creating GBM similarity matrix
    # mat = pd.read_csv("/Users/Alexis/Desktop/StuartRotation/data/simMatrix.pancan.atlas.imputed.tsv", delimiter = "\t",skiprows=0, index_col=False)
    mat = mat.set_index('sample')
    GBM_sim_matrix = pd.DataFrame(index=GBM_mRNA_ids, columns=GBM_mRNA_ids)

    #all needed comparisons
    comp=[(GBM_mRNA_ids[i], GBM_mRNA_ids[j]) for i in xrange(len(GBM_mRNA_ids)) for j in xrange(len(GBM_mRNA_ids))]
    for c in comp:
        GBM_sim_matrix.loc[c] = mat.loc[c]

    print(time.time()-startTime)

    #****takes over 4 minutes
    ### Binary Feature Creation - GBM ATRX mutations
    varianttypes = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site']
    GBM_ATRX = pd.DataFrame(0, index = range(1), columns = GBM_mRNA_ids)
    # snv_file = open("/Users/Alexis/Desktop/StuartRotation/data/snvout-2.tsv", "r")

    for line in snv_file:
        sline = line.split("\t")
        if sline[0]== 'ATRX':
            #is it a GBM sample?
            if sline[15] in GBM_snv_ids:
                #is it an impactful type?
                if sline[8] in varianttypes:
                    #add to matrix
                    GBM_ATRX[sline[15]] = 1

    print(time.time()-startTime)

    return GBM_sim_matrix, GBM_ATRX

def main():
    global startTime
    startTime = time.time()

    clinical = open("/Users/Alexis/Desktop/StuartRotation/data/clinical_out.tsv","r")
    IDmapfile = open("/Users/Alexis/Desktop/StuartRotation/data/ID_mapping.tsv", "r")
    mat = pd.read_csv("/Users/Alexis/Desktop/StuartRotation/data/simMatrix.pancan.atlas.imputed.tsv", delimiter = "\t",skiprows=0, index_col=False)
    snv_file = open("/Users/Alexis/Desktop/StuartRotation/data/snvout-2.tsv", "r")



    #Gathering GBM ATRX data
    GBM_sim_matrix, GBM_ATRX = GBMprep(clinical, IDmapfile, mat, snv_file)

    #Running Electo
    AnalyzeFeature("ATRX", GBM_sim_matrix, GBM_ATRX)

    print("done",time.time()-startTime)

if __name__ == "__main__": main()
