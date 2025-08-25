import argparse
import h3.api.numpy_int as h3
import numpy as np
import xarray as xr

try:
    input_file_default = snakemake.input[0]
    h3_resolution = snakemake.params["h3_res"]
    output_file_default = snakemake.output[0]
except NameError:
    input_file_default = None
    h3_resolution = None
    output_file_default = None

def get_h3_cells_within(ds: xr.Dataset, h3_res: int) -> np.array:
    lon_min, lon_max = ds.longitude.min().values.item(), ds.longitude.max().values.item()
    lat_min, lat_max = ds.latitude.min().values.item(), ds.latitude.max().values.item()

    lon_mid = (lon_max + lon_min) / 2

    # NOTE: h3 needs lat first.
    # NOTE: The H3 library misinterprets bounding boxes that span >180° longitude as crossing the antimeridian ("transmeridian").
    #       To avoid this, we split the longitude range into two halves if needed.
    # See: https://github.com/uber/h3/issues/1000
    # See: https://github.com/uber/h3/issues/945#issuecomment-2565448754
    bbox = [
        (lat_max, lon_min),
        (lat_max, (lon_min + lon_mid)/2),
        (lat_max, lon_mid),
        (lat_max, (lon_max + lon_mid)/2),
        (lat_max, lon_max), 
        (lat_min, lon_max),
        (lat_min, (lon_max + lon_mid)/2),
        (lat_min, lon_mid),
        (lat_min, (lon_min + lon_mid)/2),
        (lat_min, lon_min),
    ]

    poly = h3.LatLngPoly(bbox)
    return h3.h3shape_to_cells_experimental(poly, res=h3_res)

def main(input_file: str, output_file: str, h3_res: int):
    ds = xr.open_dataset(input_file, drop_variables=['time', 'surface', 'heightAboveGround'], decode_timedelta=False)

    h3_cells_int = get_h3_cells_within(ds=ds, h3_res=h3_res)
    bbox_h3_latlng = np.array([h3.cell_to_latlng(i) for i in h3_cells_int])

    coords = {"cell_ids": ("h3_cell", h3_cells_int)}

    ds_h3 = ds.interp(
        longitude=xr.DataArray(bbox_h3_latlng[:, 1], dims="h3_cell", coords=coords),
        latitude=xr.DataArray(bbox_h3_latlng[:, 0], dims="h3_cell", coords=coords),
    )
    ds_h3 = ds_h3.drop_vars(["longitude", "latitude"])
    # This is needed for xdggs
    ds_h3.cell_ids.attrs = {"grid_name": "h3", "resolution": h3_res}

    ds_h3.to_netcdf(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Aggregate ECMWF GRIB data daily.')
    parser.add_argument('--input_file', type=str, required=(input_file_default is None), default=input_file_default, help='Input GRIB file path')
    parser.add_argument('--output_file', type=str, required=(output_file_default is None), default=output_file_default, help='Output NetCDF file path')
    parser.add_argument('--h3_resolution', type=int, required=(h3_resolution is None), default=h3_resolution, help='The H3 resolution to use')

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.h3_resolution)