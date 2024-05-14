import argparse
import pandas as pd

def main(input_file):
    df = pd.read_csv(
        input_file,
        sep='\t\|\t',  
        engine='python',
        header=0,
        dtype=str,
        lineterminator='\n'
    )
    
    # Select the desired columns
    df = df[['uid', 'parent_uid', 'rank', 'sourceinfo']]

    # Extract the ncbi id from sourceinfo
    df['NCBI'] = df['sourceinfo'].str.extract(r'ncbi:(\d+)')

    df = df[['uid', 'NCBI']]

    df = df[df['NCBI'].notna() & df['NCBI'].str.strip().astype(bool)]
    df['NCBI'] = df['NCBI'].astype(int)

    print(df)

    df.to_csv('NCBI_to_ToL_file.csv', index = False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a taxonomy file.')
    parser.add_argument('input_file', type=str, help='Path to the input file')

    args = parser.parse_args()
    main(args.input_file)