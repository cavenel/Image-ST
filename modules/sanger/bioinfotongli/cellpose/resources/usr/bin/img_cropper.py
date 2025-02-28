#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

import logging
import tifffile
from aicsimageio import AICSImage
import zarr
import numpy as np

logging.basicConfig(level=logging.INFO)

def get_tile_from_tifffile(image, xmin, xmax, ymin, ymax, channel=0, zplane=0, timepoint=0, resolution_level=0):
    store = tifffile.imread(image, aszarr=True)
    zgroup = zarr.open(store, mode="r")
    if isinstance(zgroup, zarr.core.Array):
        dimension_order = zgroup.attrs["_ARRAY_DIMENSIONS"]
        if len(zgroup.shape) == 2:
            dimension_order = "YX"
        elif dimension_order == ["Q", "Q", "Q", "Y", "X"]:
            dimension_order = "QQQYX"
        else:
            logging.error(f"Unknown dimension order {image.shape}")
        image = zgroup
    else:
        image = zgroup[resolution_level]
        dimension_order = [d[0] for d in image.attrs["_ARRAY_DIMENSIONS"]]
        dimension_order = "".join(dimension_order)

    # crop = image[y_min:y_max, x_min:x_max]
    # print(image.shape)
    if dimension_order=="YX":
        tile = image[ymin:ymax, xmin:xmax]
    elif dimension_order=="YXC" or dimension_order=="YXS":
        tile = image[ymin:ymax, xmin:xmax, channel]
    elif dimension_order=="CYX" or dimension_order=="SYX":
        tile = image[channel, ymin:ymax, xmin:xmax]
    elif dimension_order=="ZYX":
        tile = image[zplane, ymin:ymax, xmin:xmax]
    elif dimension_order=="ZYXC":
        tile = image[zplane, ymin:ymax, xmin:xmax, channel]
    elif dimension_order=="YXCZ":
        tile = image[ymin:ymax, xmin:xmax, channel, zplane]
    elif dimension_order=="XYCZT":
        tile = image[ymin:ymax, xmin:xmax, channel, zplane, timepoint]
    elif dimension_order=="QQQYX":
        tile = image[0, channel, 0,  ymin:ymax, xmin:xmax]
    else:
        raise Exception(f"Unknown dimension order {dimension_order}")
    
    logging.debug(f"tile shape {tile.shape}")
    return tile

def slice_and_crop_image(image_p, x_min, x_max, y_min, y_max, zs, channels, resolution_level):
    if image_p.endswith(".tif") or image_p.endswith(".tiff"):
        crop = get_tile_from_tifffile(image_p, x_min, x_max, y_min, y_max, zplane=zs, resolution_level=resolution_level)
    else:
        # This will load the whole slice first and then crop it. So, large memroy footprint
        img = AICSImage(image_p)
        ch_ind = channels[0] if len(np.unique(channels)) == 1 else channels
        lazy_one_plane = img.get_image_dask_data(
            "ZCYX",
            T=0, # only one time point is allowed for now
            C=ch_ind,
            Z=zs)
        crop = lazy_one_plane[:, :, y_min:y_max, x_min:x_max].compute()
    return crop