#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

import logging
import tifffile
import zarr 
import numpy as np

logging.basicConfig(level=logging.INFO)

def get_tile(image, xmin, xmax, ymin, ymax, channel=0, zplane=0, timepoint=0, resolution_level=0):
    store = tifffile.imread(image, aszarr=True)
    zgroup = zarr.open(store, mode="r")
    print(zgroup) 
    if isinstance(zgroup, zarr.core.Array):
        image = np.array(zgroup)
        if len(image.shape) == 2:
            dimension_order = "YX"
        else:
            logging.error(f"Unknown dimension order {image.shape}")
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
    else:
        raise Exception(f"Unknown dimension order {dimension_order}")
    
    logging.debug(f"tile shape {tile.shape}")
    return tile