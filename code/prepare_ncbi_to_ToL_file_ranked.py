import argparse
import pandas as pd
import networkx as nx
import gc

def process_chunk(df_chunk, G, rank_cache, all_ranks):
    # Extract NCBI IDs from sourceinfo
    df_chunk['NCBI'] = df_chunk['sourceinfo'].str.extract(r'ncbi:(\d+)')
    df_chunk['uid_leaf'] = 'ott' + df_chunk['uid'].astype(str)
    df_chunk = df_chunk[df_chunk['NCBI'].notna() & df_chunk['NCBI'].str.strip().astype(bool)]
    df_chunk['NCBI'] = df_chunk['NCBI'].astype(int)
    
    # Initialize rank columns
    for rank in all_ranks:
        df_chunk[rank] = None

    total_uids = len(df_chunk['uid'])
    for i, row in df_chunk.iterrows():
        if i % 1000 == 0:
            print(f"Processing UID {i+1}/{total_uids} in chunk")
        uid = row['uid']
        
        if uid in rank_cache:
            for rank, name in rank_cache[uid].items():
                df_chunk.at[i, rank] = name
            continue
        
        ranks = {}
        try:
            ancestors = nx.ancestors(G, uid)
            for ancestor in ancestors:
                if ancestor in rank_cache:
                    ranks.update(rank_cache[ancestor])
                    break
                try:
                    parent = list(G.predecessors(ancestor))[0]
                    edge_data = G.get_edge_data(parent, ancestor)
                    rank = edge_data['rank']
                    if rank in df_chunk.columns and rank not in ranks:
                        ranks[rank] = edge_data['name']
                except (IndexError, TypeError):
                    continue
        except nx.NetworkXError:
            pass
        
        if ranks:
            for rank, name in ranks.items():
                df_chunk.at[i, rank] = name
            rank_cache[uid] = ranks
        
        if i % 10000 == 0:
            gc.collect()  # Perform garbage collection to free up memory

    df_chunk = df_chunk.drop(columns=['sourceinfo', 'uniqname', 'flags\t|'])
    return df_chunk

def main(input_file):
    print("Reading the input file...")
    chunksize = 100000  # Process the file in chunks of 100000 lines
    df_chunks = pd.read_csv(
        input_file,
        sep=r'\t\|\t',  
        engine='python',
        header=0,
        dtype=str,
        lineterminator='\n',
        chunksize=chunksize
    )

    print("Creating a directed graph...")
    G = nx.DiGraph()
    rank_cache = {}
    
    for i, df_chunk in enumerate(df_chunks):
        print(f"Processing chunk {i+1} for graph creation...")
        for _, row in df_chunk.iterrows():
            if pd.notna(row['parent_uid']):
                G.add_edge(row['parent_uid'], row['uid'], rank=row['rank'], name=row['name'])
        gc.collect()  # Perform garbage collection to free up memory
    print("Finished adding edges.")
    
    # Determine all ranks
    print("Reading the input file again to get all ranks...")
    all_ranks = set()
    for df_chunk in pd.read_csv(
        input_file,
        sep=r'\t\|\t',  
        engine='python',
        header=0,
        dtype=str,
        lineterminator='\n',
        chunksize=chunksize
    ):
        all_ranks.update(df_chunk['rank'].unique())
    all_ranks.discard('no rank')
    all_ranks.discard('no rank - terminal')
    all_ranks.discard('infraspecificname')
    all_ranks.discard('variety')
    all_ranks.discard('subphylum')
    all_ranks.discard('species subgroup')
    all_ranks.discard('subdivision')
    all_ranks.discard('superclass')
    all_ranks.discard('subspecies')
    all_ranks.discard('superphylum')
    all_ranks.discard('subtribe')
    all_ranks.discard('section')
    all_ranks.discard('species group')
    all_ranks.discard('varietas')
    all_ranks.discard('forma')
    all_ranks.discard('subkingdom')
    all_ranks.discard('subsection')
    all_ranks.discard('supertribe')
    all_ranks.discard('subterclass')
    all_ranks.discard('cohort')
    all_ranks.discard('subcohort')
    all_ranks.discard('species')
    print(f"Rank columns determined: {list(all_ranks)}")

    print("Processing and writing output in chunks...")
    output_file = 'NCBI_to_ToL_file.csv'
    with open(output_file, 'w') as f_out:
        for i, df_chunk in enumerate(pd.read_csv(
            input_file,
            sep=r'\t\|\t',  
            engine='python',
            header=0,
            dtype=str,
            lineterminator='\n',
            chunksize=chunksize
        )):
            print(f"Processing chunk {i+1} for rank information...")
            df_processed = process_chunk(df_chunk, G, rank_cache, all_ranks)
            df_processed.to_csv(f_out, index=False, header=(i == 0))
            print(f"Chunk {i+1} written to output.")
            gc.collect()  # Perform garbage collection to free up memory

    print("Finished writing all chunks to CSV file.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a taxonomy file and construct a network.')
    parser.add_argument('input_file', type=str, help='Path to the input file')

    args = parser.parse_args()
    main(args.input_file)
