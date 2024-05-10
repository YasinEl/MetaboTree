#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/code"
DATA_FOLDER = "$baseDir/data"

// param for local masst results
params.MASST_input = "$DATA_FOLDER/CCMSLIB00005465843_c27ol.csv"


// params to request masst results via usi
params.usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00006116693"
params.precursor_mz_tol = 0.05
params.mz_tol = 0.05
params.min_cos = 0.6
params.analog = false
params.analog_mass_below = 130
params.analog_mass_above = 200
params.database = "gnpsdata_index"



process DownloadData {
    publishDir "./nf_output", mode: 'copy'

    output:
    path 'redu.tsv' 

    script:
    """
    wget -q -O redu.tsv https://redu.gnps2.org/dump
    """
}

process RequestNCBIlinages {
    conda "$baseDir/envs/r_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path input_file

    output:
    path "all_redu_linages_redu.csv"

    script:
    """
    $TOOL_FOLDER/request_ncbi_linages.R $input_file all_redu_linages_redu.csv
    """
}


process CreateTaxTree {
    conda "$baseDir/envs/py_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path input_csv

    output:
    path "*.json"
    path "*.nw"

    script:
    """
    python $TOOL_FOLDER/createTxonomicTree.py --input_csv $input_csv --output_json tree_output.json --output_nw tree_output.nw
    """
}

process ProcessMASSTResults {
    conda "$baseDir/envs/r_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path redu_tsv
    path masst_csv

    output:
    path "masst_by_ncbi_output.tsv"

    script:
    """
    Rscript $TOOL_FOLDER/prepare_MASST_annotations.R $redu_tsv $masst_csv masst_by_ncbi_output.tsv
    """
}

process VisualizeTreeData {
    conda "$baseDir/envs/r_env.yml" 

    publishDir "./nf_output", mode: 'copy'

    input:
    path input_tree
    path input_lib
    path input_linage
    path input_redu

    output:
    path "tree.png"

    script:
    """
    Rscript $TOOL_FOLDER/make_ggtree.R  \
    --input_tree $input_tree \
    --input_lib $input_lib \
    --input_redu $input_redu \
    --input_lin $input_linage \
    --output_png tree.png
    """
}


process RunFastMASST {

    cache false


    input:
    val usi

    output:
    path "masst_results.json"
    path "masst_results.csv"

    script:
    """
     python $TOOL_FOLDER/fast_search.py "$usi" --precursor_mz_tol $params.precursor_mz_tol --mz_tol $params.mz_tol --min_cos $params.min_cos ${params.analog ? '--analog' : ''} --analog_mass_below $params.analog_mass_below --analog_mass_above $params.analog_mass_above --database $params.database
    """
}


workflow {
    redu = DownloadData()
    ncbi_linage = RequestNCBIlinages(redu)
    (tree_json, tree_nerwik) = CreateTaxTree(ncbi_linage)
    (masst_response_json, masst_response_csv) = RunFastMASST(params.usi)
    masst_input = ProcessMASSTResults(redu, masst_response_csv)
    VisualizeTreeData(tree_nerwik, masst_input, ncbi_linage, redu)
}