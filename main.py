from src.response_function import Radiance_Map
import os
import numpy as np
from PIL import Image as PILImage
import time
import csv
from fractions import Fraction
from pathlib import Path
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
        img = PILImage.open(path).convert("L") # grayscale for now -> böyle kaldı
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
    t_start = time.perf_counter()

    Rs, Images = examples(folder, custom, file_csv)
    t_load = time.perf_counter()
    print(f"Loading images took {t_load - t_start:.2f} s")

    radmap = Radiance_Map(Images, Rs)
    radmap.response_function()
    t_calib = time.perf_counter()
    print(f"Calibration took {t_calib - t_load:.2f} s")
    print(f"Total (load + calibration): {t_calib - t_start:.2f} s")

    output_dir = Path("outputs") / Path(folder).name
    output_dir.mkdir(parents=True, exist_ok=True)
    # plotted last, so the time you spend looking at the window isn't counted
    radmap.plot_response_function( output_dir / "response_final.png")

    # Save coefficients and exposure ratios
    radmap.save_results( output_dir / "coefficients.csv",t_calib - t_start)

    return


import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("folder")
    parser.add_argument("--custom",action="store_true",)
    parser.add_argument("--csv", default="")
    args = parser.parse_args()

    main(args.folder, args.custom, args.csv)

    # main("./images/exposure sequences/exposure sequences/virtual_camera/linear/")
    # main("./images/exposure sequences/exposure sequences/virtual_camera/noise25/")
    # main("./images/exposure sequences/exposure sequences/virtual_camera/noise50/")
    # main("./images/exposure sequences/exposure sequences/virtual_camera/noise100/")
    # main("./images/exposure sequences/exposure sequences/virtual_camera/noise150/")
    # main("./images/exposure sequences/exposure sequences/virtual_camera/srgb/")
    # main("./images/exposure sequences/exposure sequences/Canon_EOS_550D/")
    # main("./images/exposure sequences/exposure sequences/Minolta_DiMAGE_A1/")
    # main("./images/exposure sequences/exposure sequences/Nikon_Coolpix_E5400/")
    # main("./images/exposure sequences/exposure sequences/Nikon_D2H/")
    # main("./images/mine/", True, "./images/mine/exposures.csv")