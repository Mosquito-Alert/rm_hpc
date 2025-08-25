import argparse
from datetime import datetime
from ecmwf.opendata import Client

# Try to import snakemake variables, if running within Snakemake
try:
    date_default = snakemake.params["date"]
    output_file_default = snakemake.output[0]
except NameError:
    date_default = None
    output_file_default = None

def main(date: datetime, output_file: str):

    client = Client()

    client.retrieve(
        date=date,
        time=0,
        type='fc', # HRES -> Forecast
        param=['2t', 'tp'],
        step=list(range(0, 144, 3)) + list(range(144, 241, 6)),
        target=output_file
    )

if __name__ == "__main__":
    def valid_date(s):
        try:
            return datetime.strptime(s, "%Y-%m-%d").date()
        except ValueError:
            raise argparse.ArgumentTypeError(f"Not a valid date: '{s}'. Expected format: YYYY-MM-DD.")

    parser = argparse.ArgumentParser(description="Download ECMWF HRES data via CDS API.")
    parser.add_argument('--date', type=valid_date, required=(date_default is None), default=date_default, help="Date in YYYY-MM-DD format")
    parser.add_argument('--output_file', type=str, required=(output_file_default is None), default=output_file_default, help='Output filename')

    args = parser.parse_args()

    main(
        date=args.date,
        output_file=args.output_file
    )