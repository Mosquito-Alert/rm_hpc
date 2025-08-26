include: "preprocess_ecmwf_hres.smk"
include: "preprocess_gpw.smk"

rule merge_features:
    input:
        ecmwf="data/ecmwf/hres/daily_stats/{year}-{month}-{day}.nc",
        gpw="data/gpw/h3_population_density.csv",
    output:
        temp("outputs/features/{year}-{month}-{day}.nc")
    conda:
        "../envs/global.yaml"
    script:
        "../scripts/preprocess/merge_features.py"
