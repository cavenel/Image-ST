#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

"""

"""
import fire
from aicsimageio import AICSImage
from spotiflow.model import Spotiflow
import csv


def main(
        image_path:str, output_name:str, C:int,
        x_min:int, x_max:int, y_min:int, y_max:int,
        model_name:str="general",
        T:int=0,
        Z:int=0
    ):
    img = AICSImage(image_path)
    lazy_one_plane = img.get_image_dask_data(
        "YX",
        T=T,
        C=C,
        Z=Z
    )
    crop = lazy_one_plane[y_min:y_max, x_min:x_max].squeeze().compute()
    model = Spotiflow.from_pretrained(model_name)
    peaks, details  = model.predict(crop)

    with open(output_name, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['y', 'x'])  # write column names
        if len(peaks) > 0:
            # Serialize peaks to disk as CSV
            peaks[:, 0] = peaks[:, 0] + y_min
            peaks[:, 1] = peaks[:, 1] + x_min 
            writer.writerows(peaks)


if __name__ == "__main__":
    options = {
        "run" : main,
        "version" : "0.0.1"
    }
    fire.Fire(options)
