import argparse

def read_cids_from_file(file_path):
    """Read CIDs from a TSV file."""
    with open(file_path, 'r') as file:
        cids = [line.strip() for line in file]
    return cids

def create_sparql_query(pubchem_ids):
    """Generate a SPARQL query based on the provided list of PubChem IDs."""
    values_clause = " ".join([f'("{cid}")' for cid in pubchem_ids])
    query = f"""
SELECT ?chemical ?organism ?organismLabel ?ncbiTaxonomyID WHERE {{
  VALUES (?pubchemID) {{ {values_clause} }}
  ?chemical wdt:P662 ?pubchemID. # PubChem ID
  ?chemical wdt:P703 ?organism. # Found in taxon
  ?organism wdt:P685 ?ncbiTaxonomyID. # NCBI Taxonomy ID
  SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
}}
"""
    return query

def write_query_to_file(cids, output_file):
    """Write the SPARQL query to a file based on a list of PubChem IDs."""
    query = create_sparql_query(cids)
    with open(output_file, 'w') as file:
        file.write(query)
    print(f"SPARQL query for PubChem IDs has been written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate a SPARQL query for a list of PubChem IDs from a TSV file.")
    parser.add_argument("file_path", type=str, help="Path to the TSV file containing PubChem IDs.")
    parser.add_argument("cid", type=str, help="Path to the TSV file containing PubChem IDs.")
    parser.add_argument("usi", type=str, help="Path to the TSV file containing PubChem IDs.")
    args = parser.parse_args()

    file_name = 'sparql_' + args.cid + '_' + args.usi.split(":")[-1] + '.sparql'


    cids = read_cids_from_file(args.file_path)
    write_query_to_file(cids, file_name)

if __name__ == "__main__":
    main()
