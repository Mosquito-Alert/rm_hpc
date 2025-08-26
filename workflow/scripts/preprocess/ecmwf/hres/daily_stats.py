import argparse
import xarray as xr

try:
    input_file_default = snakemake.input[0]
    output_file_default = snakemake.output[0]
except NameError:
    input_file_default = None
    output_file_default = None

def main(input_file: str, output_file: str):
    ds = xr.open_dataset(input_file, decode_timedelta=False)
    ds = ds.swap_dims({'step': 'valid_time'})

    # Daily precipitation is calculated by deaccumulating the data:
    # daily_precip = accumulation_at_tomorrow_midnight - accumulation_at_today_midnight
    tp = ds.tp.sel(valid_time=ds.valid_time.dt.hour == 0).diff(dim='valid_time', label='lower')
    tp["valid_time"] = tp.valid_time.dt.date
    tp = tp.rename({"valid_time": "date"})

    # Last day won't be used due to it can not be obtaining on the deaccumulating.
    ds_valid = ds.sel(valid_time=ds.valid_time.dt.date.isin(tp.date))

    t2m_mean = ds_valid.t2m.groupby("valid_time.date").mean() - 273.15

    ds_agg = xr.Dataset(
        {
            't2m_mean': t2m_mean,
            'tp': tp
        },
    )
    ds_agg = ds_agg.assign_coords(date=ds_agg.date.astype('datetime64[s]'))

    ds_agg['t2m_mean'].attrs['title'] = 'Avg Temperature 2m'
    ds_agg['t2m_mean'].attrs['units'] = 'degrees Celsius'
    ds_agg['tp'].attrs['title'] = 'Total precipitation'
    ds_agg['tp'].attrs['units'] = 'm'  # for total precipitation

    ds_agg = ds_agg.drop_vars('step')

    ds_agg.to_netcdf(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str, required=(input_file_default is None), default=input_file_default, help='Input NetCDF file path')
    parser.add_argument('--output_file', type=str, required=(output_file_default is None), default=output_file_default, help='Output NetCDF file path')

    args = parser.parse_args()

    main(args.input_file, args.output_file)