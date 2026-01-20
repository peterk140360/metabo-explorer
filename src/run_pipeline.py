###############################################################################
# Author:      -
# Company:     -
# File:        run_pipeline.py
# created:     20.09.2025
# edited:      -
#
# Description: This script automates the full metabolite enrichment pipeline.
#              It executes the following steps sequentially:
#                1. Download and extract raw HMDB and Lipid Maps data
#                   using `collect-raw-data.py`.
#                2. Enrich metabolites with ClassyFire taxonomy
#                   using `enrich-classyfire.py`.
#                3. Further enrich metabolites with NPClassifier data
#                   using `enrich-np-classifier.py`.
#                4. Add Lipid Maps structural information
#                   using `enrich-lipidmaps.py`.
#
#              After all enrichment scripts complete successfully, the script:
#                - Dynamically locates the most recent LipidMaps-enriched JSON
#                  file in `enrichment/03LIPIDMAPS/output`.
#                - Converts the nested JSON data into a flattened pandas DataFrame,
#                  including fields `taxonomy`, `np_taxonomy`, `lm_taxonomy`, and
#                  `biological_properties`.
#                - Ensures the `lm_id` column is positioned next to `accession`.
#                - Saves the flattened DataFrame as a Parquet file in the
#                  `050_SHINY\data` folder (two levels above `src/enrichment`).
#
#              The script provides timestamped console output for each step,
#              including progress messages for running enrichment scripts and
#              JSON → Parquet conversion. The pipeline stops immediately if
#              any script fails, ensuring errors are caught early.
#
# Usage:       Run from the `src` folder:
#                  python run_pipeline.py
#
# Dependencies: Python 3.10+, pandas, numpy
#               All enrichment scripts must be present in their respective
#               folders as defined in the project structure.
#
###############################################################################

import subprocess
import sys
import os, time
from datetime import datetime
from pathlib import Path
import json
import pandas as pd
from pathlib import Path

# Root of the src folder
src_root = Path(__file__).parent.resolve()  # e.g., C:\Users\prior\Videos\050_SHINY\src

# The enrichment folder is inside src
enrichment_root = src_root / "enrichment"

# Parquet output folder two levels above enrichment (050_SHINY\data)
data_folder = enrichment_root.parent.parent / "data"
data_folder.mkdir(parents=True, exist_ok=True)


def timestamp():
    return datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")


def run_script(script_path):
    print(f"{timestamp()} Running {script_path} ...")
    try:
        result = subprocess.run([sys.executable, script_path], check=True)
        print(f"{timestamp()} Completed {script_path}\n")
    except subprocess.CalledProcessError as e:
        print(
            f"{timestamp()} ERROR: Script {script_path} failed with return code {e.returncode}"
        )
        raise e


def find_latest_json(folder, ending=".json"):
    """Return the JSON file with the most recent modification time."""
    json_files = [f for f in Path(folder).glob(f"*{ending}")]
    if not json_files:
        return None
    return max(json_files, key=lambda f: f.stat().st_mtime)


def flatten_fields(df, fields):
    for field in fields:
        if field in df.columns:
            nested_df = pd.json_normalize(df[field])
            nested_df.columns = [f"{field}_{col}" for col in nested_df.columns]
            df = pd.concat([df.drop(columns=[field]), nested_df], axis=1)
    return df


def flatten_biological_properties(df):
    if "biological_properties" in df.columns:
        bio_df = pd.json_normalize(df["biological_properties"], sep="_")
        bio_df.columns = [f"biological_properties_{col}" for col in bio_df.columns]
        df = pd.concat([df.drop(columns=["biological_properties"]), bio_df], axis=1)
    return df


if __name__ == "__main__":
    # --- 1. Run enrichment scripts sequentially ---
    scripts = [
        "collection/collect-raw-data.py",
        "enrichment/01CLASSYFIRE/enrich-classyfire.py",
        "enrichment/02NPCLASSIFIER/enrich-np-classifier.py",
        "enrichment/03LIPIDMAPS/enrich-lipidmaps.py",
    ]

    for script in scripts:
        run_script(script)

    # --- 2. Convert final LipidMaps-enriched JSON to Parquet ---
    lipidmaps_output_folder = Path("enrichment/03LIPIDMAPS/output")
    input_json = find_latest_json(lipidmaps_output_folder, ending=".json")
    if not input_json:
        raise FileNotFoundError(f"No JSON file found in {lipidmaps_output_folder}")
    print(f"{timestamp()} Found LipidMaps JSON for conversion: {input_json}")

    # Parquet output folder two levels up
    data_folder = Path(enrichment_root).parent.parent / "data"
    data_folder.mkdir(parents=True, exist_ok=True)  # Ensure the folder exists
    parquet_file = data_folder / input_json.name.replace(".json", ".parquet")
    
    # Force version suffix for testing
    # manual_suffix = "2021-11-16_hmdb_metabolites_classy_np_lm_2025-09-20"
    # parquet_file = data_folder / input_json.name.replace(".json", f"{manual_suffix}.parquet")

    # Load JSON
    with open(input_json, "r", encoding="utf-8") as f:
        data = json.load(f)

    raw_df = pd.DataFrame(data)

    # Flatten nested fields
    nested_fields = ["taxonomy", "np_taxonomy", "lm_taxonomy"]
    df = flatten_biological_properties(flatten_fields(raw_df, nested_fields))

    # Move lm_id next to accession
    if "lm_id" in df.columns and "accession" in df.columns:
        cols = list(df.columns)
        cols.remove("lm_id")
        accession_index = cols.index("accession") + 1
        cols.insert(accession_index, "lm_id")
        df = df[cols]

    # Save Parquet in 050_SHINY\data
    df.to_parquet(parquet_file, index=False)
    print(
        f"{timestamp()} Converted and flattened JSON saved as Parquet: {parquet_file}"
    )

    # Path to your app.py
    app_file = Path(__file__).parent / "app.py"

    # "Touch" the file to trigger reload
    os.utime(app_file, None)

    print(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Triggered app reload by touching {app_file}"
    )

