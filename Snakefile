rule install_deps:
    input:
        "renv.lock"
    output:
        ".deps-installed"
    shell:
        """Rscript -e 'renv::restore();file.create(".deps-installed")'"""

rule all:
    input:
        "data-plot/spag.pdf"

rule data:
    input:
        ".deps-installed",
        "data/data.R",
        "data-raw/flu raw data export from redcap 9-4-20.xlsx",
        "data-raw/HI results_HSCT_sent.xlsx",
        "data-raw/list of patients.xlsx"
    output:
        "data/data.csv"
    shell:
        "Rscript data/data.R"

rule data_plot:
    input:
        ".deps-installed",
        "data-plot/data-plot.R",
        "data/data.csv"
    output:
        "data-plot/spag.pdf"
    shell:
        "Rscript data-plot/data-plot.R"
