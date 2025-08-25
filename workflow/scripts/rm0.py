import argparse
import numpy as np
import xarray as xr

# Original code: https://github.com/mpardo1/RM_mosquito/blob/main/funcR0.R

try:
    input_file_default = snakemake.input[0]
    species_default = snakemake.params["species"]
    output_file_default = snakemake.output[0]
except NameError:
    input_file_default = None
    species_default = None
    output_file_default = None

def briere_func(scale: float, temp: float, tmin: float, tmax: float) -> float:
    # See: https://doi.org/10.1093/ee/28.1.22
    """
    Brière thermal performance curve with non-negative truncation.

    Parameters:
        temp (float): Temperature (°C).
        tmin (float): Minimum temperature (°C) for positive rate.
        tmax (float): Maximum temperature (°C) for positive rate.
        scale (float): Scaling parameter controlling maximum rate.

    Returns:
        float: Value of the Brière function, truncated at zero.
    """
    value = scale * temp * (temp - tmin) * (tmax - temp)
    return np.sqrt(value) if value > 0 else 0.0

def quad_func(scale: float, temp: float, tmin: float, tmax: float) -> float:
    """
    Brière-like quadratic function with non-negative truncation.

    Parameters:
        temp (float): Temperature (°C).
        tmin (float): Minimum temperature for positive rate.
        tmax (float): Maximum temperature for positive rate.
        scale (float): Scaling factor controlling the maximum rate.

    Returns:
        float: Value of the quadratic function, truncated at zero.
    """
    result = -scale * (temp - tmin) * (temp - tmax)
    if result < 0:
        return 0.0
    return result

def truncated_quad_func(a: float, b: float, c: float, temp: float) -> float:
    """
    Quadratic function with non-negative truncation.

    Parameters:
        temp (float): Temperature or input variable.
        a (float): Quadratic coefficient (temp^2 term).
        b (float): Linear coefficient (temp term).
        c (float): Constant term.

    Returns:
        float: Value of the quadratic function, truncated at zero.
    """
    result = a * np.power(temp, 2) + b * temp + c
    if result < 0:
        return 0.0
    return result

def hatching_func(rainfall: float, human_density: float) -> float:
    # See: https://doi.org/10.1098/rsif.2018.0761
    # Constants
    erat = 0.5
    e0 = 1.5
    evar = 0.05
    eopt = 8
    efac = 0.01
    edens = 0.01

    return (1-erat)*(
        (
            (1+e0)*np.exp(-evar*np.power(rainfall-eopt, 2))
        )/(
            np.exp(-evar*np.power(rainfall - eopt, 2)) + e0
        )
    ) + erat*(edens/(edens + np.exp(-efac*human_density)))

def rm0(species, temp: float, human_density: float, rainfall: float) -> float:
    # NOTE: in ecmwf, the total precipitation is in m, this needs mm.
    h = hatching_func(human_density=human_density, rainfall=rainfall)
    if species == 'albopictus':
        a = briere_func(scale=0.000193,tmin=10.25,tmax=38.32,temp=temp)
        f = 0.5 * briere_func(scale=0.0488,tmin=8.02,tmax=35.65,temp=temp)
        delta_adult = quad_func(scale=1.43,tmin=13.41,tmax=31.51,temp=temp)
        dE = quad_func(scale=0.00071,tmin=1.73,tmax=40.51,temp=temp)
        prob_egg2larvae = quad_func(scale=0.002663,tmin=6.668,tmax=38.92,temp=temp)
        delta_eggs = truncated_quad_func(a=0.0019328,b=-0.091868,c=1.3338874,temp=temp)
    elif species == 'aegypti':
        a = 1
        f = briere_func(scale=0.00856,tmin=14.58,tmax=34.61,temp=temp)
        delta_adult = quad_func(scale=0.004186,tmin=9.373,tmax=40.26,temp=temp)
        dE = briere_func(scale=0.0003775 ,tmin=14.88,tmax=37.42,temp=temp)
        prob_egg2larvae = quad_func(scale=0.004186,tmin=9.373,tmax=40.26,temp=temp)
        delta_eggs = truncated_quad_func(a=0.004475,b=-0.210787,c=2.552370,temp=temp)
    else:
        raise ValueError
    return np.cbrt(f * a * delta_adult * prob_egg2larvae * ((h * dE)/(h * dE + delta_eggs)))

def main(input_file, output_file, species):
    ds = xr.open_dataset(input_file)

    rm = xr.apply_ufunc(
        rm0,
        species,
        ds.t2m_mean,
        ds.population_density,
        ds.tp * 1000,  # in mm instead of m.
        input_core_dims=[],
        vectorize=True
    )

    ds_rm = xr.Dataset({"rm": rm})
    ds_rm = ds_rm.rename({'h3_cell': 'cell_ids'})
    ds_rm.cell_ids.attrs = ds.h3_cell.attrs
    ds_rm.rm.attrs['title'] = 'Mosquito basic reproductive number'
    ds_rm.to_netcdf(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str, required=(input_file_default is None), default=input_file_default, help='Input NetCDF file path')
    parser.add_argument('--output_file', type=str, required=(output_file_default is None), default=output_file_default, help='Output NetCDF file path')
    parser.add_argument('--species', type=str, required=(species_default is None), default=species_default, help='The species to use')

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.species)