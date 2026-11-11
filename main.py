from src.response_function import Radiance_Map
import os
import numpy as np
from PIL import Image as PILImage

import csv
from fractions import Fraction

def get_exposure_ratios(csv_path):
    with open(csv_path, 'r') as file:
        reader = csv.reader(file)
        row = next(reader)

    exposure_times = [float(Fraction(value.strip())) for value in row]

    R_est_list = [
        exposure_times[q] / exposure_times[q + 1]
        for q in range(len(exposure_times) - 1)
    ]

    return R_est_list

def get_images_from(folder):
    filenames = []
    for f in os.listdir(folder):
        if f.lower().endswith(('.jpg','.jpeg')):
            filenames.append(f)
    
    filenames = sorted(filenames)
    imgs = []
    for f in filenames:
        path = os.path.join(folder,f)
        img = PILImage.open(path).convert("L")
        imgs.append(np.asarray(img))
        
    if not all(img.shape == imgs[0].shape for img in imgs):
        raise ValueError("Imageların boyutları farklı")

    Image = np.stack([img.flatten() for img in imgs], axis=0)
    return Image

def examples(folder, custom, file_csv: str = ""):
    Images = get_images_from(folder)
    if not custom:
        Rs = [0.5 for _ in range(Images.shape[0] - 1)]
    else:
        Rs = get_exposure_ratios(file_csv)
    
    return Rs, Images



def main(folder, custom: bool = False, file_csv: str = ""):
    print("Radiometric Self Calibration (Mitsunaga and Nayar 1999)")
    # get the image folder path as the parameter of main
    # turn images into M, with the dimensions  P (pixel number) x Q (image count)
    Rs, Images = examples(folder, custom, file_csv)
    radmap = Radiance_Map(Images, Rs)
    radmap.response_function()
    radmap.plot_response_function()
    
    return


import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Radiometric Self Calibration (Mitsunaga and Nayar)"
    )

    parser.add_argument("folder", help="Path to the image folder")
    parser.add_argument(
        "--custom",
        action="store_true",
        help="Use custom exposure ratios from a CSV"
    )
    parser.add_argument(
        "--csv",
        default="",
        help="Path to the exposure CSV file"
    )

    args = parser.parse_args()

    main(args.folder, args.custom, args.csv)

# main("./images/exposure sequences/exposure sequences/virtual_camera/linear/")
# main("./images/exposure sequences/exposure sequences/Nikon_Coolpix_E5400/")
# main("./images/mine/", True, "./images/mine/exposures.csv")