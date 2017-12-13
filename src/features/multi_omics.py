from definitions import ROOT_DIR
from src.features.clinicaldata import ClinicalData
from src.features.genomicdata import GeneExpression, SNP, DNAMethylation, miRNAExpression, CopyNumberVariation, \
    ProteinExpression
from src.features.slide_image import WholeSlideImages


class MultiOmicsData:
    def __init__(self, cancer_type, folder_path, modalities=["WSI", "GE", "SNP", "CNV", "DNA", "MIR", "PRO"]):
        """
        Load all multi-omics TCGA data from a given folder_path with the following folder structure:

            folder_path/
                clinical/
                    genome.wustl.edu_biospecimen_sample_luad.txt
                    nationwidechildrens.org_clinical_patient_luad.txt
                gene_exp/
                    LUAD__geneExp.txt
                mirna/
                    LUAD__miRNAExp__RPM.txt
                cnv/
                    LUAD__copyNumber.txt
                protein_rppa/
                    LUAD__protein_RPPA.txt
                somatic/
                    LUAD__somaticMutation_geneLevel.txt

        :param cancer_type: TCGA cancer code name
        :param folder_path: relative directory path to the folder containing multi-omics data downloaded from TCGA-assembler
        """
        self.cancer_type = cancer_type
        self.modalities = modalities

        # LOADING DATA FROM FILES
        self.multi_omics_data = {}
        self.clinical = ClinicalData(cancer_type, folder_path + "clinical/")
        self.multi_omics_data["PATIENTS"] = self.clinical.patient
        self.multi_omics_data["BIOSPECIMENS"] = self.clinical.biospecimen

        if ("WSI" in modalities):
            self.WSI = WholeSlideImages(cancer_type, folder_path)
            self.multi_omics_data["WSI"] = self.WSI
        if ("GE" in modalities):
            self.GE = GeneExpression(cancer_type, folder_path + "gene_exp/", )
            self.multi_omics_data["GE"] = self.GE.data
        if ("SNP" in modalities):
            self.SNP = SNP(cancer_type, folder_path + "somatic/")
            self.multi_omics_data["SNP"] = self.SNP.data
        if ("MIR" in modalities):
            self.MIR = miRNAExpression(cancer_type, folder_path + "mirna/")
            self.multi_omics_data["MIR"] = self.MIR.data
        if ("DNA" in modalities):
            self.DNA = DNAMethylation(cancer_type, folder_path + "dna/")
            self.multi_omics_data["DNA"] = self.DNA.data
        if ("CNV" in modalities):
            self.CNV = CopyNumberVariation(cancer_type, folder_path + "cnv/")
            self.multi_omics_data["CNV"] = self.CNV.data
        if ("PRO" in modalities):
            self.PRO = ProteinExpression(cancer_type, folder_path + "protein_rppa/")
            self.multi_omics_data["PRO"] = self.PRO.data

        self.print_sample_sizes()

    def match_samples(self, modalities):
        """
        Return the index of bcr_sample_barcodes of the intersection of samples from all modalities

        :param modalities: An array of modalities
        :return: An pandas Index list
        """
        matched_samples = self.multi_omics_data[modalities[0]].index.copy()

        for modality in modalities:
            matched_samples = matched_samples.join(self.multi_omics_data[modality].index, how="inner")

        return matched_samples


    def normalize_data(self):
        pass

    def print_sample_sizes(self):
        for modality in self.multi_omics_data.keys():
            print(modality, self.multi_omics_data[modality].shape if hasattr(self.multi_omics_data[modality],
                                                                             'shape') else "Didn't import data")

    def get_slide_image(self):
        pass

    def get_clinical(self, ):
        pass

    def get_gene_exp(self):
        pass

    def get_miRNA_exp(self):
        pass

    def get_cnv(self):
        pass

    def get_protein_exp(self):
        pass


if __name__ == '__main__':
    folder_path = ROOT_DIR + "/data/tcga-assembler/LUAD/"
    luad_data = MultiOmicsData(cancer_type="LUAD", folder_path=folder_path,
                               modalities=["WSI", "GE", "SNP", "CNV", "MIR", "PRO"])
    print("matched samples", luad_data.match_samples(modalities=["GE", "SNP", "CNV", "MIR"]).shape)
