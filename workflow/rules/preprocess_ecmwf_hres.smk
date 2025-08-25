rule download_ecmwf_hres:
    output:
        temp("data/ecmwf/hres/raw/{year}-{month}-{day}.grib2")
    conda:
        "../envs/global.yaml"
    group:
        "ecmwf_hres_preprocess_{year}-{month}-{day}"
    params:
        date=lambda wildcards: f"{wildcards.year}-{wildcards.month}-{wildcards.day}",
    script:
        "../scripts/preprocess/ecmwf/hres/download.py"

rule convert_to_h3_ecmwf_hres:
    input:
        "data/ecmwf/hres/raw/{year}-{month}-{day}.grib2"
    output:
        temp("data/ecmwf/hres/raw_h3/{year}-{month}-{day}.nc")
    conda:
        "../envs/global.yaml"
    group:
        "ecmwf_hres_preprocess_{year}-{month}-{day}"
    resources:
        mem_mb=20480,  # 20GB
        runtime="1h"
    params:
        h3_res=config['h3_res'],
    script:
        "../scripts/preprocess/ecmwf/hres/convert_to_h3.py"

rule daily_stats_ecmwf_hres:
    input:
        "data/ecmwf/hres/raw_h3/{year}-{month}-{day}.nc"
    output:
        temp("data/ecmwf/hres/daily_stats/{year}-{month}-{day}.nc")
    conda:
        "../envs/global.yaml"
    group:
        "ecmwf_hres_preprocess_{year}-{month}-{day}"
    script:
        "../scripts/preprocess/ecmwf/hres/daily_stats.py"
