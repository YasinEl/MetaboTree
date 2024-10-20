import argparse
import pandas as pd
import os

def extract_unique_msv(file_paths):
    msv_set = set()
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split('/')
                if parts[0].startswith('MSV'):
                    msv_set.add(parts[0])
    return list(msv_set)


def USI2MASST_matchCol(x):
    elements = x.split(':')
    mzspec = elements[0]
    msvid = elements[1]
    filename = os.path.splitext(os.path.basename(elements[2]))[0]
    return f"{mzspec}:{msvid}:{filename}"

def main():
    parser = argparse.ArgumentParser(description='Process REDU data.')
    parser.add_argument('input_redu_path', type=str, help='Path to the input REDU file')
    parser.add_argument('masst1', type=str, help='Path to the input REDU file')
    parser.add_argument('masst2', type=str, help='Path to the input REDU file')
    args = parser.parse_args()
    
    # Read and process data
    dt_redu = pd.read_csv(args.input_redu_path, sep='\t')
    dt_redu = dt_redu[dt_redu['NCBITaxonomy'].str.contains('|', regex=False)]
    
    dt_redu['NCBI'] = dt_redu['NCBITaxonomy'].apply(lambda x: x.split('|')[0])
    dt_redu['tax_name'] = dt_redu['NCBITaxonomy'].apply(lambda x: x.split('|')[1])
    dt_redu = dt_redu[(dt_redu['NCBI'].notna()) & (dt_redu['NCBI'] != '') & (dt_redu['tax_name'].notna()) & (dt_redu['tax_name'] != '')]
    
    dt_redu['match_col'] = dt_redu['USI'].apply(lambda x: USI2MASST_matchCol(x))
    dt_redu = dt_redu.drop(columns=['USI'])


    #unique_msvs = extract_unique_msv([args.masst1, args.masst2])
    
    #dt_redu = dt_redu[dt_redu['ATTRIBUTE_DatasetAccession'].isin(unique_msvs)]


    dt_redu.to_csv('redu.tsv', sep='\t', index=False)

if __name__ == '__main__':
    main()
