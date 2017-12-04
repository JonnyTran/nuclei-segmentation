import numpy as np
import pandas as pd


class GenomicData:
    def __init__(self, cancer_type, file_path, columns="GeneSymbol|TCGA"):
        self.cancer_type = cancer_type

        self.data = self.preprocess_expression_table(pd.read_table(file_path, sep="\t"), columns)

        self.samples = self.data.index
        self.features = self.data.columns.tolist()
        # self.features.remove("bcr_sample_barcode")

    def preprocess_expression_table(self, df, columns):
        table = df

        # Filter columns
        table = table.filter(regex=columns)

        # Cut TCGA column names to sample barcode
        table.rename(columns=lambda x: x[:16] if ("TCGA" in x) else x, inplace=True)

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(table.columns, return_index=True)
        table = table.iloc[:, i]

        # Drop NA GeneSymbol rows
        table.dropna(axis=0, inplace=True)

        # Remove entries with unknown Gene Symbol
        table = table[table.GeneSymbol != '?']

        # Transpose dataframe to patient rows and GeneSymbol columns
        table.index = table.GeneSymbol
        table.drop(['GeneSymbol'], axis=1, inplace=True)
        table = table.T

        # Add column for patients barcode
        # table['bcr_sample_barcode'] = table.index

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(table.columns, return_index=True)
        table = table.iloc[:, i]

        return table


    def get_genes_list(self):
        return self.features

    def get_samples_list(self):
        return self.samples


class GeneExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = folder_path + "LUAD__geneExp.txt"
        super().__init__(cancer_type, file_path)


class SNP(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = folder_path + "LUAD__somaticMutation_geneLevel.txt"
        super().__init__(cancer_type, file_path)


class miRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = folder_path + "LUAD__miRNAExp__RPM.txt"
        super().__init__(cancer_type, file_path)


class CopyNumberVariation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = folder_path + "LUAD__copyNumber.txt"
        super().__init__(cancer_type, file_path)


class DNAMethylation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = folder_path + ""
        super().__init__(cancer_type, file_path)


class ProteinExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = folder_path + "LUAD__protein_RPPA.txt"
        super().__init__(cancer_type, file_path)
