import argparse
import h3
from pathlib import Path
import xarray as xr

try:
    input_file_default = snakemake.input[0]
    output_dir_default = snakemake.output[0]
except NameError:
    input_file_default = None
    output_dir_default = None

def main(input_file, output_dir):
    ds = xr.open_dataset(input_file)
    ds = ds.assign_coords(
        h3_index=xr.apply_ufunc(
            h3.int_to_str,
            ds['cell_ids'],
            vectorize=True,
            dask='parallelized',
            output_dtypes=[str]
        )
    ).set_xindex('h3_index').drop_vars('cell_ids').rename(
        {'rm': 'value'}
    )

    for date in ds.date.dt.date.values:
        df = ds.sel(date=ds.date.dt.date == date).drop_vars('date').to_dataframe()
        output_path = Path(f"{output_dir}/{date}.csv")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str, required=(input_file_default is None), default=input_file_default, help='Input NetCDF file path')
    parser.add_argument('--output_dir', type=str, required=(output_dir_default is None), default=output_dir_default, help='Output directory path')

    args = parser.parse_args()

    main(args.input_file, args.output_dir)