import pandas as pd

from src.features.clinicaldata import ClinicalData


class GenomicData:
    def __init__(self, cancer_type, file_path, clinical: ClinicalData):
        self.cancer_type = cancer_type
        self.clinical = clinical

        self.data = pd.read_table(file_path, sep="\t")
        # self.samples
        self.features = self.data.columns

    def preprocess_expression_table(self, table):
        # Cut column names to sample barcode
        table.rename(columns=lambda x: x[:16], inplace=True)

        # Drop NA rows
        table.dropna(axis=0, inplace=True)

        # Remove entries with unknown Gene Symbol
        # table = table[table.GeneSymbol != '?']

        # table.filter(regex="GeneSymbol|TCGA")

        # Get list of genes
        self.set_genes_list(list(table[table.GeneSymbol != '?']['GeneSymbol']))

        return table

    def get_genes_list(self):
        return self.genes_list

    def set_genes_list(self, genes_list):
        self.genes_list = genes_list


class GeneExpression(GenomicData):
    def __init__(self, cancer_type, folder_path, clinical: ClinicalData):
        file_path = folder_path + "LUAD__geneExp.txt"
        super().__init__(cancer_type, file_path, clinical)


class SNP(GenomicData):
    def __init__(self, cancer_type, folder_path, clinical: ClinicalData):
        file_path = folder_path + "LUAD__somaticMutation_geneLevel.txt"
        super().__init__(cancer_type, file_path, clinical)


class miRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path, clinical: ClinicalData):
        file_path = folder_path + ""
        super().__init__(cancer_type, file_path, clinical)


class CopyNumberVariation(GenomicData):
    def __init__(self, cancer_type, folder_path, clinical: ClinicalData):
        file_path = folder_path + ""
        super().__init__(cancer_type, file_path, clinical)


class DNAMethylation(GenomicData):
    def __init__(self, cancer_type, folder_path, clinical: ClinicalData):
        file_path = folder_path + ""
        super().__init__(cancer_type, file_path, clinical)


class ProteinExpression(GenomicData):
    def __init__(self, cancer_type, folder_path, clinical: ClinicalData):
        file_path = folder_path + "LUAD__protein_RPPA.txt"
        super().__init__(cancer_type, file_path, clinical)
