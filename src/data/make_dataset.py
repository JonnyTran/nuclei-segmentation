# -*- coding: utf-8 -*-
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

from definitions import ROOT_DIR


class CRCHistoPhenotypes_Patch:
    def __init__(self, patch_size=(52, 52)):
        self.classif_dir = os.path.join(ROOT_DIR, 'data/external/CRCHistoPhenotypes_2016_04_28/Classification/')
        self.patch_size = patch_size

    def load_patches(self):
        X = []
        y = []
        for img_dir in os.listdir(self.classif_dir):
            if "img" in img_dir:
                image = plt.imread(self.classif_dir + img_dir + "/" + img_dir + ".bmp")
                fibroblast = scipy.io.loadmat(self.classif_dir + img_dir + "/" + img_dir + "_fibroblast.mat")[
                    "detection"]
                epithelial = scipy.io.loadmat(self.classif_dir + img_dir + "/" + img_dir + "_epithelial.mat")[
                    "detection"]
                inflammatory = scipy.io.loadmat(self.classif_dir + img_dir + "/" + img_dir + "_inflammatory.mat")[
                    "detection"]
                others = scipy.io.loadmat(self.classif_dir + img_dir + "/" + img_dir + "_others.mat")["detection"]

                for i in range(fibroblast.shape[0]):
                    X.append(self.extract_patch(image, fibroblast[i]))
                    y.append("fibroblast")
                for i in range(epithelial.shape[0]):
                    X.append(self.extract_patch(image, epithelial[i]))
                    y.append("epithelial")
                for i in range(inflammatory.shape[0]):
                    X.append(self.extract_patch(image, inflammatory[i]))
                    y.append("inflammatory")
                for i in range(others.shape[0]):
                    X.append(self.extract_patch(image, others[i]))
                    y.append("others")
        return X, y

    def extract_patch(self, image, patch_center):
        patch_size = self.patch_size

        patch_y = int(patch_center[0] - patch_size[0] / 2)
        patch_x = int(patch_center[1] - patch_size[1] / 2)

        if patch_y < 0 or patch_x < 0 or patch_y > image.shape[0] - 1 - patch_size[0] or patch_x > image.shape[1] - 1 - \
                patch_size[1]:
            image = np.lib.pad(image, ((int(patch_size[0] / 2), int(patch_size[1] / 2)),
                                       (int(patch_size[0] / 2), int(patch_size[1] / 2)),
                                       (0, 0)), 'symmetric')
            patch_x += int(patch_size[1] / 2)
            patch_y += int(patch_size[0] / 2)

        patch_image = image[patch_x:patch_x + patch_size[0], patch_y:patch_y + patch_size[1]]
        return patch_image


if __name__ == '__main__':
    dataset = CRCHistoPhenotypes_Patch()
    X, y = dataset.load_patches()
    print(len(X), len(y))
    # pickle.dump(X, open("./CRCHistoPhenotypes_Patch_X.pickle", "wb"))
    # pickle.dump(y, open("./CRCHistoPhenotypes_Patch_y.pickle", "wb"))
