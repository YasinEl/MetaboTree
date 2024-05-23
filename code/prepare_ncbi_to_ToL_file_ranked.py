import argparse
import pandas as pd
import networkx as nx

def main(input_file):
    print("Reading the input file...")
    df = pd.read_csv(
        input_file,
        sep=r'\t\|\t',  
        engine='python',
        header=0,
        dtype=str,
        lineterminator='\n'
    )
    print("Finished reading the input file.")
    
    print("Filtering out 'no rank' and 'no rank - terminal' entries...")
    df_filtered = df[~df['rank'].isin(['no rank', 'no rank - terminal'])]
    print("Finished filtering.")
    
    print("Creating a directed graph...")
    G = nx.DiGraph()
    
    print("Adding edges to the graph...")
    for _, row in df.iterrows():
        if pd.notna(row['parent_uid']):
            G.add_edge(row['parent_uid'], row['uid'], rank=row['rank'], name=row['name'])
    print("Finished adding edges.")
    
    print("Initializing rank columns...")
    rank_columns = {rank: [] for rank in df_filtered['rank'].unique()}
    print(f"Rank columns initialized: {list(rank_columns.keys())}")
    
    print("Traversing the graph and collecting rank information...")
    total_uids = len(df['uid'])
    for i, uid in enumerate(df['uid']):
        if i % 1000 == 0:
            print(f"Processing UID {i+1}/{total_uids}...")
        ranks = {rank: None for rank in rank_columns}
        try:
            ancestors = nx.ancestors(G, uid)
            for ancestor in ancestors:
                try:
                    edge_data = G.get_edge_data(list(G.predecessors(ancestor))[0], ancestor)
                    rank = edge_data['rank']
                    if rank in ranks and ranks[rank] is None:
                        ranks[rank] = edge_data['name']
                except IndexError:
                    continue
        except nx.NetworkXError:
            pass
        for rank, name in ranks.items():
            rank_columns[rank].append(name)
    print("Finished collecting rank information.")
    
    print("Adding rank columns to the dataframe...")
    for rank, names in rank_columns.items():
        df[rank] = names
    print("Finished adding rank columns.")
    
    print("Extracting NCBI IDs from sourceinfo...")
    df['NCBI'] = df['sourceinfo'].str.extract(r'ncbi:(\d+)')
    print("Finished extracting NCBI IDs.")
    
    df['uid_leaf'] = 'ott' + df['uid'].astype(str)
    
    df = df[df['NCBI'].notna() & df['NCBI'].str.strip().astype(bool)]
    df['NCBI'] = df['NCBI'].astype(int)
    
    print("Selecting final columns...")
    final_columns = ['uid', 'NCBI', 'uid_leaf'] + list(rank_columns.keys())
    df_final = df[final_columns]
    print("Finished selecting final columns.")
    
    print("Writing output to CSV file...")
    df_final.to_csv('NCBI_to_ToL_file.csv', index=False)
    print("Finished writing to CSV file.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a taxonomy file and construct a network.')
    parser.add_argument('input_file', type=str, help='Path to the input file')

    args = parser.parse_args()
    main(args.input_file)
