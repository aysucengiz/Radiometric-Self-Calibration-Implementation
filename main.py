from src.response_function import Radiance_Map
import os
import numpy as np
from PIL import Image as PILImage

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

    Image = np.stack([img.flatten() for img in imgs], axis=1)
    return Image

def examples(folder):
    Images = get_images_from(folder)
    Rs = [2.0 for _ in range(Images.shape[1] - 1)]
    
    return Rs, Images



def main(folder):
    print("Radiometric Self Calibration (Mitsunaga and Nayar 1999)")
    # get the image folder path as the parameter of main
    # turn images into M, with the dimensions  P (pixel number) x Q (image count)
    Rs, Images = examples(folder)
    radmap = Radiance_Map(Images, Rs)
    radmap.response_function()
    radmap.plot_response_function()
    
    return