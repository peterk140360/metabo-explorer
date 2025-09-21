###############################################################################
# Author:      -
# Company:     -
# File:        enrich-classyfire.py
# Created:     20.09.2025
# Edited:      21.09.2025
#
# Description: This script parses the most recent HMDB metabolite XML file from
#              the collection folder, extracts key chemical and biological
#              information, and enriches taxonomy fields using the ClassyFire
#              API when necessary. Empty values for taxonomy or biological
#              properties are stored as None instead of empty lists.
#
#              Progress is logged in a JSON log file, and periodic backup
#              checkpoints are saved in the backup folder. The final enriched
#              dataset is written to the output folder as a JSON file, with the
#              filename derived from the original XML file and suffixed with
#              "_classy.json". For example:
#                  2021-11-17_hmdb_metabolites.xml
#              becomes:
#                  2021-11-17_hmdb_metabolites_classy.json
#
# Input:       The newest HMDB XML file in collection/download/HMDB
# Output:      Enriched JSON file in the output folder
#              Backup checkpoints in the backup folder
#              Progress log in the log folder
#
# Usage:       Run the script to parse and enrich metabolites from the latest
#              available HMDB XML dataset:
#                  python enrich-classyfire.py
#
# Notes:       - ClassyFire API is queried only when taxonomy fields are missing.
#              - Delays are included to respect API rate limits.
#              - If multiple XML files are present, the script automatically
#                selects the newest one based on modification time.
#
###############################################################################

import os
import requests
import xml.etree.ElementTree as ET
import json
import time
import datetime
import random

# ======= Config =======
# Root folder for this enrichment module
enrichment_root = os.path.dirname(os.path.abspath(__file__))

# Subfolders for logs, backups, and outputs
output_folder = os.path.join(enrichment_root, "output")
log_folder = os.path.join(enrichment_root, "log")
backup_folder = os.path.join(enrichment_root, "backup")

# Ensure folders exist
os.makedirs(output_folder, exist_ok=True)
os.makedirs(log_folder, exist_ok=True)
os.makedirs(backup_folder, exist_ok=True)

# Log file path
log_path = os.path.join(log_folder, "log.json")

# Path to the collection HMDB folder
collection_hmdb_folder = os.path.join(
    enrichment_root, "..", "..", "collection", "download", "HMDB"
)


def log_error_to_file(message):
    with open("log.txt", "a") as f:
        timestamp = datetime.datetime.now().isoformat()
        f.write(f"[{timestamp}] {message}\n")


def fetch_classyfire_taxonomy(inchikey, max_retries=5):
    url = f"http://classyfire.wishartlab.com/entities/{inchikey}.json"

    print(f"ClassyFire: querying {inchikey}")
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                time.sleep(0.1)
                return {
                    "kingdom": data.get("kingdom", {}).get("name")
                    if data.get("kingdom")
                    else None,
                    "super_class": data.get("superclass", {}).get("name")
                    if data.get("superclass")
                    else None,
                    "class": data.get("class", {}).get("name")
                    if data.get("class")
                    else None,
                    "sub_class": data.get("subclass", {}).get("name")
                    if data.get("subclass")
                    else None,
                    "direct_parent": data.get("direct_parent", {}).get("name")
                    if data.get("direct_parent")
                    else None,
                }
            elif response.status_code == 429:
                wait = min(0.5 * (2**attempt), 1.5)
                print(f"Rate limit hit for {inchikey}, retrying in {wait:.1f}s...")
                time.sleep(wait)
            else:
                print(f"HTTP {response.status_code} error for {inchikey}")
                log_error_to_file(f"HTTP error {response.status_code} for {inchikey}")
                return None

        except requests.RequestException as e:
            print(f"Request failed for {inchikey}: {e}")
            log_error_to_file(f"RequestException for {inchikey}: {str(e)}")
            time.sleep(1)

    print(f"Giving up on {inchikey} after {max_retries} attempts.")
    return None


def parse_metabolite(metabolite_element):
    ns = {"hmdb": "http://www.hmdb.ca"}
    key_order = [
        "accession",
        "status",
        "name",
        "description",
        "chemical_formula",
        "average_molecular_weight",
        "iupac_name",
        "smiles",
        "inchikey",
        "pubchem_compound_id",
        "chebi_id",
        "wikipedia_id",
    ]

    metabolite_data = {}

    # Simple flat elements — default to None (null in JSON)
    for key in key_order:
        element = metabolite_element.find(f".//hmdb:{key}", namespaces=ns)
        if element is not None and element.text and element.text.strip():
            metabolite_data[key] = element.text.strip()
        else:
            metabolite_data[key] = None

    # Taxonomy fields — always include, empty lists if missing
    taxonomy_fields = ["kingdom", "super_class", "class", "sub_class", "direct_parent"]
    taxonomy_data = {field: None for field in taxonomy_fields}

    taxonomy = metabolite_element.find("hmdb:taxonomy", namespaces=ns)
    if taxonomy is not None:
        for field in taxonomy_fields:
            el = taxonomy.find(f"hmdb:{field}", namespaces=ns)
            if el is not None and el.text and el.text.strip():
                taxonomy_data[field] = el.text.strip()

    # If taxonomy is still empty AND inchikey is valid → use ClassyFire
    if all(not v for v in taxonomy_data.values()):
        inchikey = metabolite_data.get("inchikey")
        if inchikey and inchikey != "NONE":
            classy_tax = fetch_classyfire_taxonomy(inchikey)
            if classy_tax:
                taxonomy_data.update(classy_tax)
                # Optional: sleep to be nice to the API
                time.sleep(0.1 + random.uniform(0.01, 0.03))

    metabolite_data["taxonomy"] = taxonomy_data

    # Biological properties — include tissue, cellular, biospecimen
    bio_data = {}
    bio = metabolite_element.find("hmdb:biological_properties", namespaces=ns)
    if bio is not None:
        cellular_locs = bio.find("hmdb:cellular_locations", namespaces=ns)
        biospecimen_locs = bio.find("hmdb:biospecimen_locations", namespaces=ns)
        tissue_locs = bio.find("hmdb:tissue_locations", namespaces=ns)
        # Cellular
        if cellular_locs is not None:
            cells = [
                cell.text.strip()
                for cell in cellular_locs.findall("hmdb:cellular", namespaces=ns)
                if cell.text and cell.text.strip()
            ]
            bio_data["cellular_locations"] = cells if cells else None
        else:
            bio_data["cellular_locations"] = None

        # Biospecimen
        if biospecimen_locs is not None:
            biosps = [
                biosp.text.strip()
                for biosp in biospecimen_locs.findall("hmdb:biospecimen", namespaces=ns)
                if biosp.text and biosp.text.strip()
            ]
            bio_data["biospecimen_locations"] = biosps if biosps else None
        else:
            bio_data["biospecimen_locations"] = None

        # Tissue
        if tissue_locs is not None:
            tissues = [
                tissue.text.strip()
                for tissue in tissue_locs.findall("hmdb:tissue", namespaces=ns)
                if tissue.text and tissue.text.strip()
            ]
            bio_data["tissue_locations"] = tissues if tissues else None
        else:
            bio_data["tissue_locations"] = None

    else:
        # Still include keys as empty lists
        bio_data["cellular_locations"] = None
        bio_data["biospecimen_locations"] = None
        bio_data["tissue_locations"] = None

    metabolite_data["biological_properties"] = bio_data

    return metabolite_data


def save_as_json(data, json_file):
    with open(json_file, "w") as f:
        json.dump(data, f, indent=2)


def find_input_file(folder, ending=".json"):
    files = [f for f in os.listdir(folder) if f.endswith(ending)]
    if not files:
        raise FileNotFoundError(f"No {ending} file found in {folder}")
    # Sort by modification time
    files.sort(key=lambda f: os.path.getmtime(os.path.join(folder, f)))
    return os.path.join(folder, files[-1])  # newest file


if __name__ == "__main__":
    # ======= locate input file dynamically =======
    hmdb_xml_file = find_input_file(collection_hmdb_folder, ending=".xml")
    print(f"Using latest XML file: {hmdb_xml_file}")

    # Parsing and enrichment
    metabolite_data = []
    backup_points = {1000, 10000, 50000, 100000, 150000, 200000}

    # Initialize log file
    with open(log_path, "w") as f:
        json.dump([], f)

    context = ET.iterparse(hmdb_xml_file, events=("start", "end"))
    _, root = next(context)

    for event, elem in context:
        if elem.tag.endswith("metabolite") and event == "end":
            accession = elem.findtext(
                ".//hmdb:accession",
                default="UNKNOWN",
                namespaces={"hmdb": "http://www.hmdb.ca"},
            )
            print(f"Parsing: {accession} (#{len(metabolite_data) + 1})")

            metabolite_data.append(parse_metabolite(elem))

            # Log progress every 1000
            if len(metabolite_data) % 1000 == 0:
                now = datetime.datetime.now().strftime("%H:%M:%S")
                with open(log_path, "r+") as f:
                    log = json.load(f)
                    log.append({"timestamp": now, "processed": len(metabolite_data)})
                    f.seek(0)
                    json.dump(log, f, indent=2)
                print(f"[{now}] Processed {len(metabolite_data)} metabolites...")

            # Save checkpoint
            if len(metabolite_data) in backup_points:
                backup_file = os.path.join(
                    backup_folder, f"metabolites_{len(metabolite_data)}.json"
                )
                save_as_json(metabolite_data, backup_file)
                print(f"Backup saved at {len(metabolite_data)} → {backup_file}")

            root.clear()

    # Final save using original filename + "_classy.json"
    final_file = os.path.join(
        output_folder, os.path.basename(hmdb_xml_file).replace(".xml", "_classy.json")
    )

    save_as_json(metabolite_data, final_file)
    print(f"Conversion completed. Final JSON saved at {final_file}")
