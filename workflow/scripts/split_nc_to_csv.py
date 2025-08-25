import argparse
import xarray as xr

try:
    input_file_default = snakemake.input[0]
    output_dir_default = snakemake.output[0]
except NameError:
    input_file_default = None
    output_dir_default = None

def main(input_file, output_dir):
    ds = xr.open_dataset(input_file)

    for date in ds.date.dt.date.values:
        df = ds.sel(date=ds.date.dt.date == date).rename({
            'cell_ids': 'h3_index',
            'rm': 'value'
        }).to_dataframe().reset_index()
        df.drop(columns=['date'], inplace=True)
        df.to_csv(f"{output_dir}/{date}.csv", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str, required=(input_file_default is None), default=input_file_default, help='Input NetCDF file path')
    parser.add_argument('--output_dir', type=str, required=(output_dir_default is None), default=output_dir_default, help='Output directory path')

    args = parser.parse_args()

    main(args.input_file, args.output_dir)