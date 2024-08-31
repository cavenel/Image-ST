#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

import fire
from dask_image.imread import imread
from spatialdata.models import Image2DModel
from spatialdata import SpatialData, rasterize
from spatialdata.models import ShapesModel, TableModel, PointsModel
from geopandas import GeoDataFrame
from xarray import DataArray
from spatialdata.transformations.transformations import Identity
import pandas as pd
import anndata
from shapely import from_wkt, MultiPoint, MultiPolygon
import numpy as np
from skimage.segmentation import expand_labels

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any
import logging

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)

VERSION = "0.0.1"


def load_wkts_as_shapemodel(wkt_file:str):
    with open(wkt_file, 'r') as f:
        multipoly = from_wkt(f.read())
        assert isinstance(multipoly, MultiPolygon)
    ids, polys = zip(*{i+1:poly for i, poly in enumerate(list(multipoly.geoms))}.items())
    df = pd.DataFrame({"instance_id":ids, "geometry":polys})
    df["instance_id"] = df["instance_id"].astype(int)
    df.index.name = None
    return ShapesModel.parse(GeoDataFrame(df))


def main(
        transcripts:str,
        cells_in_wkt:str,
        registered_image:str,
        out_name:str,
        pixelsize:int=1,
        y_col:str='y_location',
        x_col:str='x_location',
        feature_col:str='feature_name',
        multiscale_image:bool = True,
        raw_image_channels_to_save:list = [0, 1],
        expansion_in_pixels:int = -1,
        image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
        imread_kwargs: Mapping[str, Any] = MappingProxyType({}),):

    dapi_image = imread(registered_image, **imread_kwargs)[raw_image_channels_to_save]
    dapi_image = DataArray(dapi_image, dims=("c", "y", "x"))

    raw_image_parsed = Image2DModel.parse(
        dapi_image,
        scale_factors=[2, 2, 2, 2] if multiscale_image else None,
        transformations={"global": Identity()},
        **image_models_kwargs,
    )

    if transcripts.endswith('.csv'):
        spots = pd.read_csv(transcripts, header=0, sep=',')[[y_col, x_col, feature_col]]
    elif transcripts.endswith('.tsv'):
        spots = pd.read_csv(transcripts, header=0, sep='\t')[[y_col, x_col, feature_col]]
    elif transcripts.endswith('.wkt'):
        
        # Assuming that the wkt file contains a multipoint geometry
        with open(transcripts, 'r') as f:
            multispots = from_wkt(f.read())
        if not isinstance(multispots, MultiPoint):
            raise ValueError('Please provide a wkt file with multipoint geometry')
        spots = pd.DataFrame([(geom.y, geom.x) for geom in multispots.geoms], columns=[y_col, x_col])
        spots["feature_name"] = 'spot'
    else:
        raise ValueError('Format not recognized. Please provide a csv, tsv or wkt file')

    logger.info("Building cell shapes from wkt file")
    cell_shape = load_wkts_as_shapemodel(cells_in_wkt)

    sdata = SpatialData(
        shapes={"cell_shapes": cell_shape},
    )
    # sdata["cell_labels"] = rasterize(
    logger.info("Rasterizing cell shapes")
    cell_labels = rasterize(
        sdata["cell_shapes"],
        ["x", "y"],
        min_coordinate=[0, 0],
        max_coordinate=[dapi_image.shape[-2], dapi_image.shape[-1]],
        target_coordinate_system="global",
        target_unit_to_pixels=1.0,
    )
    sdata["dapi_image"] = raw_image_parsed

    logger.info('Assigning spots to cells')
    if expansion_in_pixels > 0:
        lab_img = expand_labels(np.array(cell_labels.data), expansion_in_pixels)
    else:
        lab_img = np.array(cell_labels.data)
    
    cell_ids = lab_img[0, 
        (spots[y_col]/pixelsize).astype(int),
        (spots[x_col]/pixelsize).astype(int)
    ]
    spots['cell_id'] = cell_ids
    logger.info("Create count matrix")
    count_matrix = spots.pivot_table(index='cell_id', columns=feature_col, aggfunc='size', fill_value=0)
    count_matrix = count_matrix.drop(count_matrix[count_matrix.index == 0].index)

    logger.info("Construct anndata object")
    adata = anndata.AnnData(X=count_matrix.values)

    REGION="cell_shapes"
    REGION_KEY = "region"
    instance_key = "instance_id"
    adata.obs[instance_key] = adata.obs.index.astype(int)
    adata.var_names = count_matrix.columns
    adata.var_names_make_unique()
    adata.obs[REGION_KEY] = REGION
    table = TableModel.parse(adata, region=REGION, region_key=REGION_KEY, instance_key=instance_key)

    points = PointsModel.parse(
        spots,
        coordinates={"x": "x_int", "y": "y_int"},
        # feature_key=instance_key,
        # instance_key=XeniumKeys.CELL_ID,
        transformations={"global": Identity()},
        # sort=True,
    )
    sdata["transcripts"] = points
    sdata["table"] = table
    sdata.write(out_name)


if __name__ == '__main__':
    options = {
        "run" : main,
        "version" : VERSION
    }
    fire.Fire(options)