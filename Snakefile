rule install_deps:
    input:
        "renv.lock"
    output:
        ".deps-installed"
    shell:
        """Rscript -e 'renv::restore();file.create(".deps-installed")'"""

rule all:
    input:
        "report/report.pdf",
        "pm-flu-rct.zip"

rule data:
    input:
        ".deps-installed",
        "data/data.R",
        "data-raw/flu raw data export from redcap 9-4-20.xlsx",
        "data-raw/HI results_HSCT_sent.xlsx",
        "data-raw/list of patients.xlsx"
    output:
        "data/data.csv",
        "data/ili.csv",
        "data/seroprotection.csv",
        "data/seroprotection_combined.csv",
        "data/seroconversion.csv",
        "data/seroconversion_combined.csv",
        "data/titre.csv",
        "data/adverse_events.csv"
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
        "data/data.csv",
        "data/adverse_events.csv"
    output:
        "data-table/nobs.tex",
        "data-table/nobs.csv",
        "data-table/prop-ili.csv",
        "data-table/prop-ili.tex",
        "data-table/prop-seroprotection.csv",
        "data-table/prop-seroprotection.tex",
        "data-table/prop-seroconversion.csv",
        "data-table/prop-seroconversion.tex",
        "data-table/prop-seroprotection_combined.csv",
        "data-table/prop-seroprotection_combined.tex",
        "data-table/prop-seroconversion_combined.csv",
        "data-table/prop-seroconversion_combined.tex",
        "data-table/mid-est.tex",
        "data-table/mid-est.csv",
        "data-table/mid-est-pvals.tex",
        "data-table/mid-est-pvals.csv",
        "data-table/gmr-pvals.csv",
        "data-table/gmr-pvals.tex",
        "data-table/gmr.csv",
        "data-table/gmr.tex",
        "data-table/adverse_events.tex",
        "data-table/adverse_events.csv",
        "data-table/adverse_events_summary.tex",
        "data-table/adverse_events_summary.csv",
        "data-table/prop-seroprotection-pvals.tex",
        "data-table/prop-seroprotection-pvals.csv",
        "data-table/prop-seroconversion-pvals.tex",
        "data-table/prop-seroconversion-pvals.csv",
        "data-table/prop-seroprotection_combined-pvals.tex",
        "data-table/prop-seroprotection_combined-pvals.csv",
        "data-table/prop-seroconversion_combined-pvals.tex",
        "data-table/prop-seroconversion_combined-pvals.csv"
    shell:
        "Rscript data-table/data-table.R"

rule fit:
    input:
        ".deps-installed",
        "data/read_data.R",
        "fit/fit.R",
        "data/data.csv"
    output:
        "fit/fits-ili.csv",
        "fit/fits-titre.csv",
        "fit/fits-seroprotection.csv",
        "fit/fits-seroprotection_combined.csv",
        "fit/fits-seroconversion.csv",
        "fit/fits-seroconversion_combined.csv",
        "fit/fit-interpret-ili.tex",
        "fit/fit-interpret-titre.tex",
        "fit/fit-interpret-seroprotection.tex",
        "fit/fit-interpret-seroprotection_combined.tex",
        "fit/fit-interpret-seroconversion.tex",
        "fit/fit-interpret-seroconversion_combined.tex",
        "fit/formula-ili.tex",
        "fit/formula-titre.tex",
        "fit/formula-seroprotection.tex",
        "fit/formula-seroprotection_combined.tex",
        "fit/formula-seroconversion.tex",
        "fit/formula-seroconversion_combined.tex",
        "fit/vars-ili.tex",
        "fit/vars-titre.tex",
        "fit/vars-seroprotection.tex",
        "fit/vars-seroprotection_combined.tex",
        "fit/vars-seroconversion.tex",
        "fit/vars-seroconversion_combined.tex",
        "fit/ili-interpretation.csv",
        "fit/titre-interpretation.csv",
        "fit/seroprotection-interpretation.csv",
        "fit/seroprotection_combined-interpretation.csv",
        "fit/seroconversion-interpretation.csv",
        "fit/seroconversion_combined-interpretation.csv"
    shell:
        "Rscript fit/fit.R"

rule fit_table:
    input:
        ".deps-installed",
        "fit-table/fit-table.R",
        "fit/fits-ili.csv",
        "fit/fits-titre.csv"
    output:
        "fit-table/fit-table-titre.tex",
        "fit-table/fit-table-titre.csv",
        "fit-table/fit-table-ili.tex",
        "fit-table/fit-table-ili.csv",
        "fit-table/fit-table-seroprotection.tex",
        "fit-table/fit-table-seroprotection.csv",
        "fit-table/fit-table-seroprotection_combined.tex",
        "fit-table/fit-table-seroprotection_combined.csv",
        "fit-table/fit-table-seroconversion.tex",
        "fit-table/fit-table-seroconversion.csv",
        "fit-table/fit-table-seroconversion_combined.tex",
        "fit-table/fit-table-seroconversion_combined.csv",
        "fit-table/fit-table-titre-pval.tex",
        "fit-table/fit-table-titre-pval.csv",
        "fit-table/fit-table-ili-pval.tex",
        "fit-table/fit-table-ili-pval.csv",
        "fit-table/fit-table-seroprotection-pval.tex",
        "fit-table/fit-table-seroprotection-pval.csv",
        "fit-table/fit-table-seroprotection_combined-pval.tex",
        "fit-table/fit-table-seroprotection_combined-pval.csv",
        "fit-table/fit-table-seroconversion-pval.tex",
        "fit-table/fit-table-seroconversion-pval.csv",
        "fit-table/fit-table-seroconversion_combined-pval.tex",
        "fit-table/fit-table-seroconversion_combined-pval.csv"
    shell:
        "Rscript fit-table/fit-table.R"

rule report:
    input:
        ".deps-installed",
        "report/report.tex",
        "data-table/nobs.tex",
        "data-table/mid-est.tex",
        "data-table/prop-ili.tex",
        "fit-table/fit-table-titre.tex",
        "fit-table/fit-table-ili.tex",
        "fit/fit-interpret-titre.tex",
        "fit/fit-interpret-ili.tex",
        "fit/formula-titre.tex",
        "fit/formula-ili.tex",
        "fit/vars-titre.tex",
        "fit/vars-ili.tex",
        "data-plot/spag.pdf"
    output:
        "report/report.pdf",
        temp("report/report.log"),
        temp("report/report.synctex.gz"),
        temp("report/report.fls"),
        temp("report/report.fdb_latexmk"),
        temp("report/report.aux"),
        temp("report/report.bbl"),
        temp("report/report.blg"),
        temp("report/report.out")
    shell:
        "cd report && latexmk -synctex=1 -interaction=nonstopmode -file-line-error -pdf report.tex"

rule zip:
    input:
        "report/report.pdf"
    output:
        "pm-flu-rct.zip"
    shell:
        "zip -r pm-flu-rct.zip . -x 'renv/library*' '.snakemake*' '.deps-installed'"
