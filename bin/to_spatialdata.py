#!/usr/bin/env python3
import fire
from dask_image.imread import imread
from spatialdata.models import Image2DModel
from spatialdata import SpatialData, rasterize
from spatialdata.models import ShapesModel, TableModel
from geopandas import GeoDataFrame
from shapely import from_wkt
from xarray import DataArray
from spatialdata.transformations.transformations import Identity

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any

VERSION = "0.0.1"


def main(
        label_image:str,
        count_matrix:str,
        wkt:str,
        out_name:str,
        multiscale_image:bool = True,
        registered_image:str=None,
        image_models_kwargs: Mapping[str, Any] = MappingProxyType({}),
        imread_kwargs: Mapping[str, Any] = MappingProxyType({}),):

    if registered_image:
        raw_image = imread(registered_image, **imread_kwargs)
        raw_image = DataArray(raw_image, dims=("c", "y", "x"))


        raw_image_parsed = Image2DModel.parse(
            raw_image,
            scale_factors=[2, 2, 2, 2] if multiscale_image else None,
            transformations={"global": Identity()},
            **image_models_kwargs,
        )

    with open (wkt, "r") as fh:
        shapes=from_wkt(fh.read())

    ids, geoms = [], []
    for i, geom in enumerate(shapes.geoms):
        ids.append(i+1)
        geoms.append(geom)

    df = pd.DataFrame({'instance_id': ids, 'geometry': geoms})
    tab = GeoDataFrame(df)
    shape = ShapesModel.parse(tab)
    sdata = SpatialData(
        shapes={"shapes": shape},
    )
    sdata["label"] = rasterize(
        sdata["shapes"],
        ["x", "y"],
        min_coordinate=[0, 0],
        max_coordinate=[20399, 20399],
        target_coordinate_system="global",
        target_unit_to_pixels=1.0,
    )

    SpatialData(
        images={"raw_image": raw_image_parsed},
    ).write(out_name)


if __name__ == '__main__':
    options = {
        "run" : main,
        "version" : VERSION
    }
    fire.Fire(options)