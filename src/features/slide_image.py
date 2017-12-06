import os

import large_image


class WholeSlideImages:
    def __init__(self, cancer_type, dir_path):
        self.cancer_type = cancer_type

        for file in os.listdir(dir_path)[0:2]:
            if file.endswith(".svs"):
                print(file)
                imagePath = os.path.join(dir_path, file)
                source = large_image.getTileSource(imagePath)
                thumbnail, mimeType = source.getThumbnail(
                    width=1024, height=1024, encoding='JPEG')
                open('./thumbnail.jpg', 'wb').write(thumbnail)


if __name__ == '__main__':
    wsi = WholeSlideImages("LUAD", "/media/jonny_admin/540GB/Research/TCGA_LUAD-WSI/")
