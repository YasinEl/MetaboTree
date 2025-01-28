import argparse
import pandas as pd
from pathlib import Path

def generate_sparql(smiles_list, group_name, search_mode):
    smiles_values = ' '.join(f'"{smiles}"' for smiles in smiles_list)
    query = f"""#title: Which organisms contain indolic scaffolds? Count occurrences, group and order the results by the parent taxon.
PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX p: <http://www.wikidata.org/prop/>
PREFIX idsm: <https://idsm.elixir-czech.cz/sparql/endpoint/>
SELECT DISTINCT ?parent_taxon ?parent_taxon_name  ?taxon ?taxon_name ?ncbi ?compound ?inchikey  WHERE {{
  SERVICE idsm:wikidata {{
    VALUES ?SUBSTRUCTURE {{
      {smiles_values} # smiles scaffolds
    }}
    ?compound sachem:substructureSearch [
      sachem:query ?SUBSTRUCTURE;
      sachem:searchMode sachem:{search_mode};
      sachem:chargeMode sachem:defaultChargeAsAny;
      sachem:isotopeMode sachem:defaultIsotopeAsStandard;
      sachem:aromaticityMode sachem:aromaticityDetectIfMissing;
      sachem:stereoMode sachem:ignoreStereo;
      sachem:tautomerMode sachem:ignoreTautomers;
      sachem:radicalMode sachem:ignoreSpinMultiplicity;
      sachem:topn '-1'^^xsd:integer;
      sachem:internalMatchingLimit '1000000'^^xsd:integer
    ].
  }}
  hint:Prior hint:runFirst "true"^^xsd:boolean.
  ?compound p:P703 ?statement;
    wdt:P235 ?inchikey.
  ?statement ps:P703 ?taxon.
  ?taxon wdt:P225 ?taxon_name;
    wdt:P171 ?parent_taxon;
    wdt:P685 ?ncbi.
  ?parent_taxon wdt:P225 ?parent_taxon_name.
  SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
}}
"""
    return query

def main():
    parser = argparse.ArgumentParser(description="Generate SPARQL queries for each group in a TSV file.")
    parser.add_argument("input_tsv", help="Path to the input TSV file containing SMILES and grouping variable.")
    parser.add_argument("output_dir", help="Directory to save the generated SPARQL files.")
    
    args = parser.parse_args()
    input_tsv = Path(args.input_tsv)
    output_dir = Path(args.output_dir)

    if not input_tsv.exists():
        raise FileNotFoundError(f"Input TSV file '{input_tsv}' does not exist.")

    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    # Load the TSV file
    df = pd.read_csv(input_tsv, sep="\t")


    if "smiles" not in df.columns or "name" not in df.columns or "match" not in df.columns:
        raise ValueError("The TSV file must contain 'smiles', 'name', and 'match' columns.")

    # Group by the grouping variable and process each group
    for group_name, group_data in df.groupby("name"):
        unique_smiles = group_data["smiles"].dropna().unique()
        if len(unique_smiles) == 0:
            continue

        # Map 'match' values to the appropriate search mode
        match_value = group_data["match"].iloc[0]
        search_mode = "substructureSearch" if match_value == "substructure" else "exactSearch"

        sparql_query = generate_sparql(unique_smiles, group_name, search_mode)
        output_file = output_dir / f"sparql_{group_name}.sparql"

        # Save the SPARQL query to a file
        with open(output_file, "w") as f:
            f.write(sparql_query)

        print(f"Generated SPARQL file for name '{group_name}' at '{output_file}'")

if __name__ == "__main__":
    main()
