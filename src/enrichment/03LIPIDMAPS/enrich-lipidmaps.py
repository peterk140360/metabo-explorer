###############################################################################
# Author:      -
# Company:     -
# File:        enrich-lipidmaps.py
# created:     20.09.2025
# edited:      21.09.2025
#
# Description: This script enriches HMDB metabolites with structural and
#              classification information from the LIPID MAPS Structure
#              Database (LMSD). It combines information from an SDF file
#              containing lipid structures with metabolites previously
#              enriched with NPClassifier data.
#
# Workflow:
#   1. Dynamically loads the newest JSON output from the NPClassifier enrichment
#      (located in 02NPCLASSIFIER/output).
#   2. Dynamically loads the LMSD SDF file from collection/download/LIPIDMAPS.
#   3. Parses the SDF to index lipids by INCHI_KEY and other identifiers
#      (PUBCHEM_CID, CHEBI_ID, NAME, FORMULA, SYSTEMATIC_NAME, SMILES).
#   4. Enriches each metabolite entry:
#       - First attempts a direct match using INCHI_KEY.
#       - If no direct match, performs fallback matching using combinations
#         of secondary identifiers (requires at least 2 matching criteria).
#       - Tracks ambiguous matches as errors.
#       - Populates `lm_id` and `lm_taxonomy` fields for each metabolite.
#   5. Logs each match, fallback, and no-match situation in a log file.
#   6. Saves the enriched metabolite dataset to the output folder in JSON format.
#
# Input:
#   - NPClassifier JSON: enrichment/02NPCLASSIFIER/output/<dynamic>.json
#   - LMSD SDF file: collection/download/LIPIDMAPS/<dynamic>.sdf
#
# Output:
#   - Enriched JSON file: enrichment/03LIPIDMAPS/output/<input>_lm_<version>.json
#   - Log file: enrichment/03LIPIDMAPS/log/log.txt
#
# Notes:
#   - Dynamic detection of input JSON and SDF makes the script robust to
#     filename changes.
#   - Handles missing or empty fields by storing `None`.
#   - Ambiguous fallback matches are logged for review.
###############################################################################

import os
import json
from rdkit import Chem
from datetime import datetime
from collections import Counter

# ======= Config =======
# Root folder for this enrichment module
enrichment_root = os.path.dirname(os.path.abspath(__file__))

# Paths for NP classifier output and .sdf file
npclassifier_output_folder = os.path.join(enrichment_root, "..", "02NPCLASSIFIER", "output")
lipidmaps_folder = os.path.join(enrichment_root, "..", "..", "collection", "download", "LIPIDMAPS")

# Output and log folders
output_folder = os.path.join(enrichment_root, "output")
log_folder = os.path.join(enrichment_root, "log")
os.makedirs(output_folder, exist_ok=True)
os.makedirs(log_folder, exist_ok=True)


def timestamp():
    return datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")


def log(msg, logfile=None):
    print(msg)
    if logfile:
        with open(logfile, "a", encoding="utf-8") as f:
            f.write(msg + "\n")


def find_input_file(folder, ending=".json"):
    files = [f for f in os.listdir(folder) if f.endswith(ending)]
    if not files:
        raise FileNotFoundError(f"No {ending} file found in {folder}")
    # Sort by modification time
    files.sort(key=lambda f: os.path.getmtime(os.path.join(folder, f)))
    return os.path.join(folder, files[-1])  # newest file


if __name__ == "__main__":
    # ======= Config =======
    enrichment_root = os.path.dirname(os.path.abspath(__file__))

    # Input JSON from NPClassifier
    input_json = find_input_file(npclassifier_output_folder, ending=".json")
    print(f"Found NPClassifier JSON file: {input_json}")

    # Input SDF file from LIPID MAPS
    sdf_file = find_input_file(lipidmaps_folder, ending=".sdf")
    print(f"Using latest SDF file: {sdf_file}")

    # Extract version from SDF filename (before "_structures.sdf")
    sdf_basename = os.path.basename(sdf_file)
    sdf_version = sdf_basename.replace("_structures.sdf", "")
    print(f"Detected LIPID MAPS version: {sdf_version}")

    # Output and log folders
    output_folder = os.path.join(enrichment_root, "output")
    log_folder = os.path.join(enrichment_root, "log")
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(log_folder, exist_ok=True)

    # Output file paths -> append _lm_{version}.json
    output_json = os.path.join(
        output_folder,
        os.path.basename(input_json).replace(".json", f"_lm_{sdf_version}.json")
    )
    log_file = os.path.join(log_folder, "log.txt")

    # --- Init log ---
    with open(log_file, "w", encoding="utf-8") as f:
        f.write(f"{timestamp()} Starting enrichment process...\n")

    # --- Load SDF into indexed lipid_data ---
    log(f"{timestamp()} Loading SDF file: {sdf_file}", log_file)
    lipid_data = {}
    index_maps = {
        "pubchem": {},
        "chebi": {},
        "name": {},
        "formula": {},
        "systematic_name": {},
        "smiles": {}
    }
    suppl = Chem.SDMolSupplier(sdf_file)

    for mol in suppl:
        if mol is None:
            continue

        inchikey = mol.GetProp('INCHI_KEY') if mol.HasProp("INCHI_KEY") else None
        if not inchikey:
            continue

        lipid_entry = {
            "lm_id": mol.GetProp('LM_ID') if mol.HasProp("LM_ID") else None,
            "category": mol.GetProp('CATEGORY') if mol.HasProp("CATEGORY") else None,
            "main_class": mol.GetProp('MAIN_CLASS') if mol.HasProp("MAIN_CLASS") else None,
            "sub_class": mol.GetProp('SUB_CLASS') if mol.HasProp("SUB_CLASS") else None,
            "smiles": mol.GetProp('SMILES') if mol.HasProp("SMILES") else None,
            "pubchem": mol.GetProp('PUBCHEM_CID') if mol.HasProp("PUBCHEM_CID") else None,
            "chebi": mol.GetProp('CHEBI_ID') if mol.HasProp("CHEBI_ID") else None,
            "name": mol.GetProp('NAME') if mol.HasProp("NAME") else None,
            "formula": mol.GetProp('FORMULA') if mol.HasProp("FORMULA") else None,
            "systematic_name": mol.GetProp('SYSTEMATIC_NAME') if mol.HasProp("SYSTEMATIC_NAME") else None
        }

        lipid_data[inchikey] = lipid_entry

        for key in index_maps:
            val = lipid_entry[key]
            if val:
                index_maps[key].setdefault(val, []).append(inchikey)

    log(f"{timestamp()} Parsed {len(lipid_data)} lipid entries from SDF.", log_file)

    # --- Load JSON ---
    log(f"{timestamp()} Loading JSON file: {input_json}", log_file)
    with open(input_json, "r", encoding="utf-8") as f:
        metabolites = json.load(f)
    log(f"{timestamp()} Loaded {len(metabolites)} metabolite entries.", log_file)

    # --- Enrich each entry ---
    matches = 0
    misses = 0
    errors = 0
    fallback_counters = {k: 0 for k in ["INCHI_KEY","PUBCHEM_CID","NAME","CHEBI_ID","FORMULA","SYSTEMATIC_NAME","SMILES"]}
    fallback_variants = Counter()

    for idx, entry in enumerate(metabolites):
        inchi_key = entry.get("inchikey")
        smiles = entry.get("smiles")
        pubchem = entry.get("pubchem_cid")
        chebi = entry.get("chebi_id")
        name = entry.get("name")
        formula = entry.get("chemical_formula")
        systematic_name = entry.get("iupac_name")

        match = None
        reason = "no match"

        if inchi_key and inchi_key in lipid_data:
            match = lipid_data[inchi_key]
            reason = "matched by INCHI_KEY"
            fallback_counters["INCHI_KEY"] += 1
            fallback_variants["INCHI_KEY"] += 1
        else:
            criteria = []
            possible_keys = set()

            if pubchem and pubchem in index_maps["pubchem"]:
                possible_keys.update(index_maps["pubchem"][pubchem])
                criteria.append("PUBCHEM_CID")
            if name and name in index_maps["name"]:
                possible_keys.update(index_maps["name"][name])
                criteria.append("NAME")
            if chebi and chebi in index_maps["chebi"]:
                possible_keys.update(index_maps["chebi"][chebi])
                criteria.append("CHEBI_ID")
            if formula and formula in index_maps["formula"]:
                possible_keys.update(index_maps["formula"][formula])
                criteria.append("FORMULA")
            if systematic_name and systematic_name in index_maps["systematic_name"]:
                possible_keys.update(index_maps["systematic_name"][systematic_name])
                criteria.append("SYSTEMATIC_NAME")
            if smiles and smiles in index_maps["smiles"]:
                possible_keys.update(index_maps["smiles"][smiles])
                criteria.append("SMILES")

            if len(criteria) >= 2 and len(possible_keys) == 1:
                match_key = list(possible_keys)[0]
                match = lipid_data[match_key]
                reason = f"fallback match by: {', '.join(criteria)}"
                for c in criteria:
                    fallback_counters[c] += 1
                fallback_variants[reason] += 1
            elif len(criteria) >= 2 and len(possible_keys) > 1:
                errors += 1
                log(f"{timestamp()} ERROR: Ambiguous fallback match for {entry.get('accession', f'IDX_{idx}')} using {', '.join(criteria)}", log_file)

        if match:
            matches += 1
            entry["lm_id"] = match["lm_id"]
            entry["lm_taxonomy"] = {
                "category": match["category"],
                "main_class": match["main_class"],
                "sub_class": match["sub_class"]
            }
            log(f"{timestamp()} MATCH ({reason}): {entry.get('accession', f'IDX_{idx}')}", log_file)
        else:
            misses += 1
            entry["lm_id"] = None
            entry["lm_taxonomy"] = {
                "category": None,
                "main_class": None,
                "sub_class": None
            }
            log(f"{timestamp()} NO MATCH: {entry.get('accession', f'IDX_{idx}')}", log_file)

    # --- Save JSON ---
    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(metabolites, f, indent=2, ensure_ascii=False)

    # --- Summary ---
    log(f"{timestamp()} Enrichment complete.", log_file)
    log(f"  Total entries: {len(metabolites)}", log_file)
    log(f"  Matches:       {matches}", log_file)
    log(f"  No matches:    {misses}", log_file)
    log(f"  Errors (ambiguous fallback): {errors}", log_file)
    log(f"  Output:        {output_json}", log_file)
    log("   Match breakdown:", log_file)
    for key, count in fallback_counters.items():
        log(f"    {key}: {count}", log_file)
    log("  Fallback variant breakdown:", log_file)
    for variant, count in fallback_variants.items():
        log(f"    {count} fallback match by: {variant}", log_file)
