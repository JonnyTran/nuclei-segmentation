from definitions import ROOT_DIR
from src.features.clinicaldata import ClinicalData
from src.features.genomicdata import GeneExpression, SNP, DNAMethylation, miRNAExpression, CopyNumberVariation, \
    ProteinExpression


class MultiOmicsData:
    def __init__(self, cancer_type, folder_path):
        """
        Load all multi-omics TCGA data from given :param folder_path: with the following folder structure:


        :param cancer_type:
        :param folder_path:
        """
        self.cancer_type = cancer_type
        self.clinical = ClinicalData(cancer_type, folder_path + "clinical/")
        # self.slide_images = TissueSlideImages(cancer_type, folder_path, clinical=self.clinical)
        self.gene_exp = GeneExpression(cancer_type, folder_path + "gene_exp/", clinical=self.clinical)
        self.snp = SNP(cancer_type, folder_path + "somatic/", clinical=self.clinical)
        self.miRNA_exp = miRNAExpression(cancer_type, folder_path + "mirna/", clinical=self.clinical)
        self.dna_methyl = DNAMethylation(cancer_type, folder_path + "dna/", clinical=self.clinical)
        self.cnv = CopyNumberVariation(cancer_type, folder_path + "cnv/", clinical=self.clinical)
        self.protein_exp = ProteinExpression(cancer_type, folder_path + "protein/", clinical=self.clinical)

    def match_data(self):
        pass

    def normalize_data(self):
        pass

    def get_slide_image(self):
        pass

    def get_clinical(self):
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
    # table = pd.read_table(ROOT_DIR+"/data/tcga-assembler/LUAD/clinical/nationwidechildrens.org_clinical_patient_luad.txt", sep="\t")
    folder_path = "/data/tcga-assembler/LUAD/"
    luad_data = MultiOmicsData(cancer_type="LUAD", folder_path=ROOT_DIR + folder_path)
