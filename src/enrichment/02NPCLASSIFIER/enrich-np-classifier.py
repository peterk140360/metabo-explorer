###############################################################################
# Author:      -
# Company:     -
# File:        enrich-np-classifier.py
# created:     20.09.2025
# edited:      22.09.2025
#
# Description: This script enriches HMDB metabolite data with natural product
#              classification information using the NPClassifier API.
#              It automatically loads the newest JSON output produced by the
#              ClassyFire enrichment (located in 01CLASSYFIRE/output) and
#              processes each metabolite entry in parallel using multiple
#              threads.
#
#              For each metabolite:
#                - The SMILES string is extracted and safely URL-encoded.
#                - If the SMILES is valid, the NPClassifier API is queried
#                  using the encoded string.
#                - NP taxonomy fields ("pathway", "super_class", "class") are
#                  populated based on the API results.
#                - Empty or missing fields are consistently stored as `None`.
#                - Failed queries or invalid SMILES are logged and tracked.
#
#              Progress is logged in a text file and error types are counted.
#              The final enriched dataset is written to the "output" folder,
#              and a detailed log is saved in the "log" folder.
#
# Input:       Newest JSON file from ClassyFire enrichment located in:
#              enrichment/01CLASSYFIRE/output
#
# Output:      Enriched JSON file in:
#              enrichment/02NPCLASSIFIER/output
#              Log file in:
#              enrichment/02NPCLASSIFIER/log
#
# Usage:       Run the script to classify metabolites using NPClassifier.
#              The script dynamically detects the ClassyFire JSON file and
#              appends "_np" to generate the output filename.
#
# Notes:       - Parallel processing is used to speed up API queries.
#              - MAX_WORKERS can be adjusted for the number of threads.
#              - The unwrap function ensures empty lists are converted to None.
#              - All SMILES strings are URL-encoded before being sent to the API
#                to prevent server errors due to reserved characters.
###############################################################################

import os
import json
import requests
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from collections import Counter
from time import time

# ======= Config =======
# Root folder for this enrichment module
enrichment_root = os.path.dirname(os.path.abspath(__file__))

# Paths for ClassyFire output (input for NP classifier)
classyfire_output_folder = os.path.join(enrichment_root, "..", "01CLASSYFIRE", "output")

# Output and log folders for NP classifier
output_folder = os.path.join(enrichment_root, "output")
log_folder = os.path.join(enrichment_root, "log")
os.makedirs(output_folder, exist_ok=True)
os.makedirs(log_folder, exist_ok=True)

# Maximum number of threads for parallel processing
MAX_WORKERS = 20


def unwrap(value):
    if isinstance(value, list):
        return value[0] if value else None
    return value


def timestamp():
    return datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")


def classify(index, entry):
    if not isinstance(entry, dict):
        return (
            index,
            {"np_taxonomy": {"pathway": None, "super_class": None, "class": None}},
            f"{timestamp()} Index {index}: entry is not a dict",
            False,
        )

    smiles = entry.get("smiles")
    accession = entry.get("accession", f"INDEX_{index}")

    # If SMILES is missing or not a string, write empty np_taxonomy
    if not isinstance(smiles, str) or not smiles.strip():
        entry["np_taxonomy"] = {"pathway": None, "super_class": None, "class": None}
        return (
            index,
            entry,
            f"{timestamp()} {accession}: No valid SMILES (null or empty)",
            False,
        )

    smiles = smiles.strip()
    encoded_smiles = urllib.parse.quote(smiles, safe="") # Encoded SMILES for URL

    try:
        url = f"https://npclassifier.gnps2.org/classify?smiles={encoded_smiles}"
        response = requests.get(url, timeout=30)

        if response.status_code == 200:
            result = response.json()
            entry["np_taxonomy"] = {
                "pathway": unwrap(result.get("pathway_results")),
                "super_class": unwrap(result.get("superclass_results")),
                "class": unwrap(result.get("class_results")),
            }
            return index, entry, "", True
        else:
            entry["np_taxonomy"] = {"pathway": None, "super_class": None, "class": None}
            return (
                index,
                entry,
                f"{timestamp()} {accession}: HTTP error {response.status_code}",
                False,
            )

    except Exception as e:
        entry["np_taxonomy"] = {"pathway": None, "super_class": None, "class": None}
        return index, entry, f"{timestamp()} {accession}: Exception - {str(e)}", False


def find_input_file(folder, ending=".json"):
    files = [f for f in os.listdir(folder) if f.endswith(ending)]
    if not files:
        raise FileNotFoundError(f"No {ending} file found in {folder}")
    # Sort by modification time
    files.sort(key=lambda f: os.path.getmtime(os.path.join(folder, f)))
    return os.path.join(folder, files[-1])  # newest file


if __name__ == "__main__":
    # ======= locate input file dynamically =======
    input_file = find_input_file(classyfire_output_folder, ending=".json")
    print(f"Found ClassyFire JSON file: {input_file}")

    # Output and log files
    output_file = os.path.join(
        output_folder, os.path.basename(input_file).replace(".json", "_np.json")
    )
    log_file = os.path.join(log_folder, "log.txt")

    # ======= start timer =======
    start_time = time()
    error_counter = Counter()

    # ======= load data =======
    with open(input_file, "r", encoding="utf-8") as f:
        metabolites = json.load(f)
    print(f"{timestamp()} Loaded {len(metabolites)} entries from {input_file}")

    # ======= ordered parallel processing =======
    results = [None] * len(metabolites)
    log_lines = []
    success_count = 0
    fail_count = 0

    print(
        f"{timestamp()} Starting parallel classification with {MAX_WORKERS} threads...\n"
    )

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(classify, i, entry): i
            for i, entry in enumerate(metabolites)
        }
        for i, future in enumerate(as_completed(futures), start=1):
            index, entry, log_line, success = future.result()
            results[index] = entry
            print(
                f"{timestamp()} Result received: {entry.get('accession')} - {'Success' if success else 'Failed'}"
            )

            if log_line:
                log_lines.append(log_line)
                # Count unique error types
                parts = log_line.split(": ", 2)
                if len(parts) == 3:
                    error_msg = parts[2]
                    error_counter[error_msg] += 1

            if success:
                success_count += 1
            else:
                fail_count += 1

    # ======= save results =======
    with open(output_file, "w", encoding="utf-8") as f_out:
        json.dump(results, f_out, indent=2, ensure_ascii=False)

    with open(log_file, "w", encoding="utf-8") as f_log:
        for line in log_lines:
            f_log.write(line + "\n")

    # ======= summary =======
    elapsed = time() - start_time
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)

    print(f"\n{timestamp()} Processing complete.")
    print(f"  Time taken: {minutes} min {seconds} sec")
    print(f"  Total entries: {len(metabolites)}")
    print(f"  Successful: {success_count}")
    print(f"  Failed:     {fail_count}")
    print(f"  Output written to: {output_file}")
    print(f"  Log saved to:      {log_file}")

    if error_counter:
        print("\n  Error breakdown:")
        for msg, count in error_counter.most_common():
            print(f"    - {msg}: {count}")
