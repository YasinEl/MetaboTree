import argparse
import requests
import logging
import json
import pandas as pd
import time

def _fast_masst(params):
    URL = "https://fasst.gnps2.org/search"
    for attempt in range(5):
        try:
            search_api_response = requests.post(URL, data=params, timeout=300)
            search_api_response.raise_for_status()
            logging.debug("fastMASST response={}".format(search_api_response.status_code))
            return search_api_response.json()
        except requests.exceptions.RequestException as e:
            logging.error(f"Request failed on attempt {attempt + 1}: {str(e)}")
            if attempt < 4:  # 4 because attempt starts from 0 and we need 5 attempts
                time.sleep(1)
            else:
                return {}

def fast_masst(params):
    try:
        response_json = _fast_masst(params)
        if not response_json:
            return pd.DataFrame(columns=['Delta Mass', 'USI', 'Charge', 'Cosine', 'Matching Peaks', 'Dataset'])
        return response_json
    except Exception as e:
        logging.error(f"Error processing the MASST request: {str(e)}")
        return pd.DataFrame(columns=['Delta Mass', 'USI', 'Charge', 'Cosine', 'Matching Peaks', 'Dataset'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fetch MASST search results from GNPS")
    parser.add_argument("usi_or_lib_id", help="USI or Library ID for the search")
    parser.add_argument("--precursor_mz_tol", type=float, default=0.02, help="Precursor M/Z tolerance")
    parser.add_argument("--mz_tol", type=float, default=0.02, help="Fragment M/Z tolerance")
    parser.add_argument("--min_cos", type=float, default=0.7, help="Minimum cosine similarity")
    parser.add_argument("--analog", action='store_true', help="Enable analog search")
    parser.add_argument("--analog_mass_below", type=int, default=130, help="Analog mass below")
    parser.add_argument("--analog_mass_above", type=int, default=200, help="Analog mass above")

    args = parser.parse_args()

    databases = ["gnpsdata_index", "gnpsdata_index_11_25_23"]
    all_results = []

    for db in databases:
        params = {
            "usi": args.usi_or_lib_id,
            "library": db,
            "analog": "Yes" if args.analog else "No",
            "delta_mass_below": args.analog_mass_below,
            "delta_mass_above": args.analog_mass_above,
            "pm_tolerance": args.precursor_mz_tol,
            "fragment_tolerance": args.mz_tol,
            "cosine_threshold": args.min_cos,
        }

        result = fast_masst(params)
        if isinstance(result, pd.DataFrame) and not result.empty:
            result['Database'] = db
            all_results.append(result)

    if all_results:
        # Combining results from all databases
        combined_df = pd.concat(all_results, ignore_index=True)
        # Drop duplicates
        combined_df.drop_duplicates(subset=['USI'], inplace=True)

        # Ensure the expected columns are present
        expected_columns = ['Delta Mass', 'USI', 'Charge', 'Cosine', 'Matching Peaks', 'Dataset']
        for col in expected_columns:
            if col not in combined_df.columns:
                combined_df[col] = pd.NA

        # Selecting specific columns if needed
        combined_df = combined_df[expected_columns]
    else:
        # Create an empty DataFrame with the required columns if no results were fetched
        combined_df = pd.DataFrame(columns=['Delta Mass', 'USI', 'Charge', 'Cosine', 'Matching Peaks', 'Dataset'])

    # Saving the DataFrame to a CSV file
    combined_df.to_csv('masst_results.csv', index=False)

    print("CSV file has been created successfully.")
