import argparse
import pandas as pd
import xarray as xr

try:
    ecmwf_file_default = snakemake.input['ecmwf']
    gpw_file_default = snakemake.input['gpw']
    output_file_default = snakemake.output[0]
except NameError:
    ecmwf_file_default = None
    gpw_file_default = None
    output_file_default = None

def main(ecmwf_file, gpw_file, output_file):
    ds_ecmwf = xr.open_dataset(ecmwf_file).set_index(h3_cell="cell_ids")
    ds_gpw = xr.Dataset.from_dataframe(
        pd.read_csv(gpw_file, index_col='cell')
    ).rename({
        'cell': 'h3_cell',
    })

    ds_gpw = ds_gpw.assign_coords(
        h3_cell=ds_gpw.h3_cell.astype("uint64")
    )

    ds_combined = xr.merge([ds_ecmwf, ds_gpw], join='inner')
    ds_combined.to_netcdf(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate GPW to H3.")
    parser.add_argument("--ecmwf_file", type=str, required=(ecmwf_file_default is None), default=ecmwf_file_default, help="Path to ECMWF input file")
    parser.add_argument("--gpw_file", type=str, required=(gpw_file_default is None), default=gpw_file_default, help="Path to GPW input file")
    parser.add_argument("--output_file", type=str, required=(output_file_default is None), default=output_file_default, help="Path to output file")

    args = parser.parse_args()
    main(
        ecmwf_file=args.ecmwf_file,
        gpw_file=args.gpw_file,
        output_file=args.output_file
    )