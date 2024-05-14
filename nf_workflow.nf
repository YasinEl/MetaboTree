#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/code"
DATA_FOLDER = "$baseDir/data"


// keys
params.email = "$baseDir/keys/email.txt"
params.entrez_key = "$baseDir/keys/entrez_key.txt"


// param for local masst results
params.MASST_input = "$DATA_FOLDER/CCMSLIB00000204966_capsacin.csv"


// params for tree type 
params.tree_type = 'treeoflife' //taxonomic

// params for SarQL
params.pubchemid = "6140"

// params to request masst results via usi
params.usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00006116693"
params.precursor_mz_tol = 0.05
params.mz_tol = 0.05
params.min_cos = 0.6
params.analog = false
params.analog_mass_below = 130
params.analog_mass_above = 200



process DownloadData {
    conda "$baseDir/envs/py_env.yml"

    publishDir "./nf_output", mode: 'copy'

    output:
    path 'redu.tsv' 

    script:
    """
    wget -q -O redu.tsv https://redu.gnps2.org/dump
    python $TOOL_FOLDER/prepare_redu_table.py redu.tsv
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

process Process_MASST_Wikidata_Pubchem_Results {
    conda "$baseDir/envs/py_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path redu_tsv
    path masst_csv
    path sparql_csv
    path ncbi_response_csv

    output:
    path "masst_by_ncbi_output.tsv"
    path "sparql_by_ncbi_output.tsv"
    path "ncbi_by_ncbi_output.tsv"

    script:
    """
    python $TOOL_FOLDER/process_masst_wikidata_pubchem_results.py $redu_tsv $masst_csv $sparql_csv $ncbi_response_csv
    """
}

process RunFastMASST {
    conda "$baseDir/envs/py_env.yml"
    
    cache false


    input:
    val usi

    output:
    //path "masst_results.json"
    path "masst_results.csv"

    script:
    """
    python $TOOL_FOLDER/fast_search.py "$usi" --precursor_mz_tol $params.precursor_mz_tol --mz_tol $params.mz_tol --min_cos $params.min_cos ${params.analog ? '--analog' : ''} --analog_mass_below $params.analog_mass_below --analog_mass_above $params.analog_mass_above 
    """
}


process RunWikidataSparql {
    

    input:
    path pubchemids

    publishDir "./nf_output", mode: 'copy'

    output:
    path "sparql_query.csv"

    script:
    """
     python $TOOL_FOLDER/make_sparql.py  $pubchemids
     wikidata-dl sparql_query.sparql
     mv wikidata/sparql_query.csv sparql_query.csv
    """
}


process MakeTreeRings {
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

process getNCBIRecords {
    cache false

    conda "$baseDir/envs/py_env.yml"


    publishDir "./nf_output", mode: 'copy'

    input:
    path pubchemids

    output:
    path "ncbi_records.csv"

    script:
    """
    python $TOOL_FOLDER/get_ncbi_records.py $DATA_FOLDER/biosystems_taxonomy $DATA_FOLDER/biosystems_pcsubstance $pubchemids
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
    path input_sparql
    path input_ncbi

    output:
    path "check_this.csv"
    path "tree.png"

    script:
    """
    Rscript $TOOL_FOLDER/make_ggtree.R  \
    --input_tree $input_tree \
    --input_lib $input_lib \
    --input_redu $input_redu \
    --input_lin $input_linage \
    --input_sparql $input_sparql \
    --input_ncbi $input_ncbi \
    --output_png tree.png
    """
}


process GetStereoCIDs {
    conda "$baseDir/envs/py_env.yml"

    input:
    val pubchemid

    output:
    path "stereoisomers.tsv"

    script:
    """
    echo "Getting stereo pubchem IDs"
    python $TOOL_FOLDER/get_stereo_cids.py  $pubchemid
    """
}


process Request_treeoflife_ids {

    cache false

    conda "$baseDir/envs/py_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path input_csv
    path ncbi_to_ToL_file

    output:
    path "modified_tree_of_life.csv"

    script:
    """
    python $TOOL_FOLDER/request_redu_from_tree_of_life.py $input_csv $ncbi_to_ToL_file
    """
}


process Update_ReDU_df_with_TOL {
    conda "$baseDir/envs/py_env.yml"


    publishDir "./nf_output", mode: 'copy'

    input:
    path treeoflife_response
    path redu

    output:
    path "redu_updated.csv"

    script:
    """
    python $TOOL_FOLDER/add_tree_of_life_info_to_redu.py --json_path $treeoflife_response --redu_path $redu
    """
}


process CreateTOLTree {
    conda "$baseDir/envs/py_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path input_csv

    output:
    path "tree.nw"

    script:
    """
    python $TOOL_FOLDER/create_tree_of_life.py $input_csv $DATA_FOLDER/labelled_supertree.tre
    """
}

process Make_ID_table {
    conda "$baseDir/envs/py_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path redu
    path sparql
    path ncbi

    output:
    path "all_organism_table.csv"

    script:
    """
    python $TOOL_FOLDER/extend_organsims_from_databases.py $redu $sparql $ncbi
    """
}



process prepare_ncbi_to_ToL_file {

    conda "$baseDir/envs/py_env.yml"

    output:
    path "NCBI_to_ToL_file.csv"

    script:
    """
    python $TOOL_FOLDER/prepare_ncbi_to_ToL_file.py $DATA_FOLDER/taxonomy.tsv
    """
}


workflow {



    cids = GetStereoCIDs(params.pubchemid)
    sparql_response_csv = RunWikidataSparql(cids)
    ncbi_response_csv = getNCBIRecords(cids)


    redu = DownloadData()
    redu_w_datasets = Make_ID_table(redu, sparql_response_csv, ncbi_response_csv)



    if (params.tree_type == 'taxonomic'){
        ncbi_linage = RequestNCBIlinages(redu_w_datasets)
        (tree_json, tree_nerwik) = CreateTaxTree(ncbi_linage)
    }
    if (params.tree_type == 'treeoflife'){
        // format NCBI -> ToL file
        ncbi_to_ToL_file = prepare_ncbi_to_ToL_file()

        redu_w_datasets = Request_treeoflife_ids(redu_w_datasets, ncbi_to_ToL_file)
        // redu_w_datasets = Update_ReDU_df_with_TOL(ncbi_to_tol, redu_w_datasets)
        tree_nerwik = CreateTOLTree(redu_w_datasets)
    }



    (masst_response_csv) = RunFastMASST(params.usi)
    //(kingdom, superclass) = MakeTreeRings(redu)






    // (masst_results, sparql_results, ncbi_results) = Process_MASST_Wikidata_Pubchem_Results(redu, masst_response_csv, sparql_response_csv, ncbi_response_csv)
//    (masst_results, sparql_results, ncbi_results) = Process_MASST_Wikidata_Pubchem_Results(redu, params.MASST_input, sparql_response_csv, ncbi_response_csv)
//    VisualizeTreeData(tree_nerwik, masst_results, ncbi_linage, redu, sparql_results, ncbi_results)
}