import argparse
import h3
from h3ronpy.pandas.raster import raster_to_dataframe
import rasterio

try:
    input_file_default = snakemake.input[0]
    h3_resolution = snakemake.params["h3_res"]
    output_file_default = snakemake.output[0]
except NameError:
    input_file_default = None
    h3_resolution = None
    output_file_default = None

def main(input_file: str, output_file: str, h3_res: int):
    with rasterio.open(input_file) as src:
        df = raster_to_dataframe(
            src.read(1),
            src.transform,
            h3_resolution=h3_res,
            nodata_value=src.nodata,
            compact=False,
            geo=False,
        )
        # Rename 'value' column
        df = df.rename(columns={'value': 'population_density'})
        # Export to CSV
        df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate GPW to H3.")
    parser.add_argument("--input_file", type=str, required=(input_file_default is None), default=input_file_default, help="Path to input TIF")
    parser.add_argument('--h3_resolution', type=int, required=(h3_resolution is None), default=h3_resolution, help='The H3 resolution to use')
    parser.add_argument("--output_file", type=str, required=(output_file_default is None), default=output_file_default, help="Path to output CSV file")

    args = parser.parse_args()
    main(
        input_file=args.input_file,
        output_file=args.output_file,
        h3_res=args.h3_resolution
    )