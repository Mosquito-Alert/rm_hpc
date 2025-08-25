rule convert_to_h3_gpw:
    input:
        "data/gpw/gpw_v4_population_density_rev11_2020_2pt5_min.tif",
    output:
        "data/gpw/h3_population_density.csv"
    conda:
        "../envs/global.yaml"
    params:
        h3_res=config['h3_res'],
    script:
        "../scripts/preprocess/gpw/convert_to_h3.py"

