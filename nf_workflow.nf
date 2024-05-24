#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/code"
DATA_FOLDER = "$baseDir/data"



// download https://pypi.org/project/itolapi/1.0/#files for tree vis 

// keys
params.email = "$baseDir/keys/email.txt"
params.entrez_key = "$baseDir/keys/entrez_key.txt"


// param for local masst results
params.MASST_input = "$DATA_FOLDER/CCMSLIB00000204966_capsacin.csv"

// param 
params.file = ''

// params for tree type 
params.tree_type = 'treeoflife' //taxonomic

// params for SarQL
params.pubchemid = "x"
params.molname = "chlorophyll"

// params to request masst results via usi
params.usi = "$DATA_FOLDER/test_usis.tsv"
params.precursor_mz_tol = 0.05
params.mz_tol = 0.05
params.min_cos = 0.6
params.analog = false
params.analog_mass_below = 130
params.analog_mass_above = 200



process PrepareReDU {
    conda "$baseDir/envs/py_env.yml"

    output:
    path 'redu.tsv' 

    script:
    """
    python $TOOL_FOLDER/prepare_redu_table.py "$DATA_FOLDER/redu.tsv" "$DATA_FOLDER/masst_datatsets_filenames.txt" "$DATA_FOLDER/masst_datatsets_filenames2.txt"
    """
}

process RequestNCBIlinages {
    conda "$baseDir/envs/r_env.yml"

    input:
    path input_file

    output:
    path "all_redu_linages_redu.csv"

    script:
    """
    $TOOL_FOLDER/request_ncbi_linages.R $input_file "$DATA_FOLDER/all_redu_linages.csv" all_redu_linages_redu.csv
    """
}


process CreateTaxTree {
    conda "$baseDir/envs/py_env.yml"

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

    publishDir "./nf_output", mode: 'copy'
    
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

    conda "$baseDir/envs/py_env.yml"
    

    input:
    each pubchemids

    publishDir "./nf_output", mode: 'copy'

    output:
    path "sparql_*.csv"

    script:
    """
     python $TOOL_FOLDER/make_sparql.py  $pubchemids $params.pubchemid $params.molname
     wikidata-dl *.sparql
     mv wikidata/*.csv ./
    """
}


process MakeTreeRings {
    conda "$baseDir/envs/r_env.yml"

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
    conda "$baseDir/envs/py_env.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    each pubchemids

    output:
    path "ncbi_records.csv"

    script:
    """
    python $TOOL_FOLDER/get_ncbi_records.py $DATA_FOLDER/biosystems_taxonomy $DATA_FOLDER/biosystems_pcsubstance $pubchemids
    """
}



process VisualizeTreeData {
    cache false

    conda "$baseDir/envs/r_env.yml" 

    publishDir "./nf_output", mode: 'copy'

    input:
    each input_tree
    each input_masst
    each input_redu
    each mol_plot
    each mol_name

    output:
    path "*.png"
    path "*.tsv", optional: true

    script:
    """
    Rscript $TOOL_FOLDER/make_ggtree_plot.R  \
    --input_tree $input_tree \
    --input_masst $input_masst \
    --input_redu $input_redu \
    --mol_plot $mol_plot \
    --mol_name $mol_name \
    --usi $params.usi \
    --cid $params.pubchemid \
    --output_png tree.png
    """
}


process GetStereoCIDs {
    
    conda "$baseDir/envs/py_env.yml"
    
    input:
    each pubchemid

    output:
    path "stereoisomers.tsv"
    path "main_name.tsv"

    script:
    """
    python $TOOL_FOLDER/get_stereo_cids.py $pubchemid
    """
}


process plot_molecule {
    
    conda "$baseDir/envs/py_env.yml"
    
    input:
    each pubchemid

    output:
    path "molecule.png"

    script:
    """
    python $TOOL_FOLDER/plot_structure.py  $pubchemid 
    """
}



process Request_treeoflife_ids {


    conda "$baseDir/envs/py_env.yml"

    input:
    each input_redu
    each input_linage
    each ncbi_to_ToL_file

    output:
    path "modified_tree_of_life.csv"
    path "linage_table.csv"

    script:
    """
    python $TOOL_FOLDER/request_redu_from_tree_of_life.py $input_redu $input_linage $ncbi_to_ToL_file
    """
}


process Update_ReDU_df_with_TOL {
    conda "$baseDir/envs/py_env.yml"


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
    each input_csv

    output:
    path "*.nw"

    script:
    """
    python $TOOL_FOLDER/create_tree_of_life.py --input_csv $input_csv --tree_path $DATA_FOLDER/labelled_supertree.tre --usi $params.usi --cid $params.pubchemid 
    """
}

process Make_ID_table {

    conda "$baseDir/envs/py_env.yml"

    input:
    path redu
    val sparql
    // each ncbi

    output:
    path "all_organism_table.csv"

    script:
    """
    python $TOOL_FOLDER/extend_organsims_from_databases.py $redu $sparql 
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


process generateEmpressPlot {

    // conda "$baseDir/envs/qiime2-amplicon-2024.2-py38-linux-conda.yml"

    conda "$baseDir/envs/env_empress_custom.yml"

    publishDir "./nf_output", mode: 'copy'

    input:
    path tree
    path masst_results

    output:
    path "**"

    script:
    """
    empress tree-plot --tree  $tree   --feature-metadata  $masst_results  --output-dir empress_results
    """
}



process makeTreeAnnotationTable {
    cache false

    conda "$baseDir/envs/r_env.yml" 

    publishDir "./nf_output", mode: 'copy'

    input:
    each input_tree
    each input_masst
    each input_redu
    each input_linage

    output:
    path "*.tsv"

    script:
    """
    Rscript $TOOL_FOLDER/make_tree_annotations.R  \
    --input_tree $input_tree \
    --input_masst $input_masst \
    --input_redu $input_redu \
    --input_linage $input_linage \
    --usi $params.usi \
    --cid $params.pubchemid \
    """
}




workflow run_phyloMASST {

    (cids, mol_name) = GetStereoCIDs(params.pubchemid)
    mol_plot = plot_molecule(params.pubchemid)
    sparql_response_csv = RunWikidataSparql(cids)
    // ncbi_response_csv = getNCBIRecords(cids)


    redu = PrepareReDU()
    redu_w_datasets = Make_ID_table(redu, sparql_response_csv)//, ncbi_response_csv)

    ncbi_linage = RequestNCBIlinages(redu_w_datasets)

    if (params.tree_type == 'taxonomic'){
        // ncbi_linage = RequestNCBIlinages(redu_w_datasets)
        (tree_json, tree_nerwik) = CreateTaxTree(ncbi_linage)
    }
    if (params.tree_type == 'treeoflife'){
        ncbi_to_ToL_file = prepare_ncbi_to_ToL_file() // format NCBI -> ToL file

        (redu_w_datasets, ncbi_linage) = Request_treeoflife_ids(redu_w_datasets, ncbi_linage, ncbi_to_ToL_file)
        tree_nerwik = CreateTOLTree(redu_w_datasets)
    }
    masst_response_csv = RunFastMASST(params.usi)
    //(kingdom, superclass) = MakeTreeRings(redu)



    (ggtree_plot, masst_results) = VisualizeTreeData(tree_nerwik, masst_response_csv, redu_w_datasets, mol_plot, mol_name)  

    

    makeTreeAnnotationTable(tree_nerwik, masst_response_csv, redu_w_datasets, ncbi_linage)

    // generateEmpressPlot(tree_nerwik, masst_results)

    
}


workflow {
    run_phyloMASST()
}