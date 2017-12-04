import pandas as pd

from definitions import ROOT_DIR


class ClinicalData:
    clinical_patient_colsname = ['bcr_patient_barcode', 'gender', 'race', 'histologic_diagnosis.1',
                                 'ajcc_pathologic_tumor_stage',
                                 ]
    biospecimen_sample_colsname = ['bcr_sample_barcode', 'sample_type']

    def __init__(self, cancer_type, folder_path):
        self.cancer_type = cancer_type
        self.patient = pd.read_table(folder_path + "nationwidechildrens.org_clinical_patient_luad.txt",
                                     sep="\t",
                                     skiprows=[1, 2],
                                     na_values="[Not Available]",
                                     usecols=ClinicalData.clinical_patient_colsname
                                     )

        self.biospecimen = pd.read_table(folder_path + "genome.wustl.edu_biospecimen_sample_luad.txt",
                                         sep="\t",
                                         skiprows=[1, ],
                                         na_values="[Not Available]",
                                         usecols=ClinicalData.biospecimen_sample_colsname
                                         )

        self.patient_barcodes = self.patient["bcr_patient_barcode"].tolist()
        self.sample_barcodes = self.biospecimen["bcr_sample_barcode"].tolist()

    def get_patient_barcodes(self):
        return self.patient_barcodes

    def get_sample_barcodes(self):
        return self.sample_barcodes


if __name__ == '__main__':
    # table = pd.read_table(ROOT_DIR+"/data/tcga-assembler/LUAD/clinical/nationwidechildrens.org_clinical_patient_luad.txt", sep="\t")
    folder_path = "/data/tcga-assembler/LUAD/clinical/"
    luad_clinical = ClinicalData(cancer_type="LUAD", folder_path=ROOT_DIR + folder_path)
