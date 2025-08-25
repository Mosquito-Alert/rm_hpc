# RM0 - Pipeline HPC

This repository contains the code for running a high-performance computing (HPC) pipeline that get the mosquito basic reproductive number based on https://github.com/mpardo1/RM_mosquito

The pipeline is orchestrated using [Snakemake](https://snakemake.github.io) which automates tasks related to data preparation, model training, and prediction.

## Data sources
The model needs data from the following data sources:
- Realtime [ECMWF - HRES](https://confluence.ecmwf.int/display/DAC/ECMWF+open+data%3A+real-time+forecasts+from+IFS+and+AIFS)
- [GPWv4 - Population Density](https://doi.org/10.7927/H49C6VHW)

## Data preparation

### Download GPWv4 - Population Density

Please download the file `gpw-v4-population-density-rev11_2020_2pt5_min_tif.zip` export it and place it in the following path: `data/gpw/gpw_v4_population_density_rev11_2020_2pt5_min.tif`

## Deployment

### 1.Install Snakemake

We recommend using [Mamba](https://mamba.readthedocs.io/en/latest/) (a faster drop-in replacement for Conda). If you don't have Conda or Mamba installed, consider installing [Miniforge](https://github.com/conda-forge/miniforge).

Install Snakemake, Snakedeploy, and necessary plugins:

```bash
mamba create -c conda-forge -c bioconda --name snakemake snakemake=9.8.1 snakedeploy=0.11.0
```

If you're running on an HPC with SLURM, install additional plugins:
```bash
mamba install -n snakemake -c bioconda snakemake-executor-plugin-slurm=1.5.0
```

Activate the environment:

```bash
conda activate snakemake
```

### 2. Deploy the workflow

Create and move into a project directory:
```bash
mkdir -p path/to/project-workdir
cd path/to/project-workdir
```

Deploy the workflow using Snakedeploy:
```bash
snakedeploy deploy-workflow <URL_TO_THIS_REPO> . --tag <DESIRED_TAG>
```

This will create two directories:
- `workflow/`: contains the deployed Snakemake module
- `config/`: contains configuration files

### 3. Configure workflow

Edit `config/config.yaml` to specify your settings (paths, parameters, etc.) according to your data and environment

### 4. Run workflow

#### Local execution with conda

```bash
snakemake --cores all --sdm conda
```

#### HPC execution with SLURM
Use the provided SLURM profile:
```bash
snakemake --cores all --sdm conda --profile slurm
```

For advanced features such as cluster execution, cloud deployments, and workflow customization, see the [Snakemake documentation](https://snakemake.readthedocs.io/).