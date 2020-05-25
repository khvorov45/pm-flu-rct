rule install_deps:
    input:
        "renv.lock"
    output:
        ".deps-installed"
    shell:
        """Rscript -e 'renv::restore();file.create(".deps-installed")'"""

rule all:
    input:
        "report/pm-flu-rct.zip"

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
        "data/read_data.R",
        "data-plot/data-plot.R",
        "data/data.csv"
    output:
        "data-plot/spag.pdf"
    shell:
        "Rscript data-plot/data-plot.R"

rule data_table:
    input:
        ".deps-installed",
        "data/read_data.R",
        "data-table/data-table.R",
        "data/data.csv"
    output:
        "data-table/nobs.tex"
    shell:
        "Rscript data-table/data-table.R"

rule fit:
    input:
        ".deps-installed",
        "data/read_data.R",
        "fit/fit.R",
        "data/data.csv"
    output:
        "fit/fits.csv"
    shell:
        "Rscript fit/fit.R"

rule fit_table:
    input:
        ".deps-installed",
        "fit-table/fit-table.R",
        "fit/fits.csv"
    output:
        "fit-table/fit-table.tex",
        "fit-table/fit-interpret.tex"
    shell:
        "Rscript fit-table/fit-table.R"

rule report:
    input:
        ".deps-installed",
        "report/report.tex",
        "data-table/nobs.tex",
        "fit-table/fit-table.tex",
        "fit-table/fit-interpret.tex",
        "data-plot/spag.pdf"
    output:
        "report/report.pdf",
        temp("report/report.log"),
        temp("report/report.synctex.gz"),
        temp("report/report.fls"),
        temp("report/report.fdb_latexmk"),
        temp("report/report.aux")
    shell:
        "cd report && latexmk -synctex=1 -interaction=nonstopmode -file-line-error -pdf report.tex"

rule zip:
    input:
        "report/report.pdf"
    output:
        "report/pm-flu-rct.zip"
    shell:
        "zip -r report/pm-flu-rct.zip . -x renv/library/**\* .git .snakemake"
