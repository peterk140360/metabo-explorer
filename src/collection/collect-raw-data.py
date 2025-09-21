###############################################################################
# Author:      -
# Company:     -
# File:        collect-raw-data.py
# created:     19.09.2025 18:36
# edited:      20.09.2025
#
# Description: The script downloads the HMDB metabolites dataset (ZIP file)
#              and the LIPID MAPS Structure Database (LMSD ZIP file), storing
#              them in separate subfolders under "collection/download" (HMDB and
#              LIPIDMAPS). It then extracts both ZIP files into their respective
#              folders. Extracted files are renamed to include the remote version
#              date scraped from the official download pages.
#
# Example:     2021-11-17_hmdb_metabolites.xml
#              2025-09-20_structures.sdf
#
# Resources:   https://hmdb.ca/downloads
#              https://lipidmaps.org/databases/lmsd/download
#
# To-Do:       Update url in line 43 to /current/hmdb_metabolites.zip
#
###############################################################################

import os
import sys
import requests
from tqdm import tqdm
import zipfile
import shutil
import urllib.request
import re

# ======= Config =======
# Root folders
collection_folder = os.path.dirname(os.path.abspath(__file__))
download_folder = os.path.join(collection_folder, "download")

# Subfolders for datasets
hmdb_folder = os.path.join(download_folder, "HMDB")
lm_folder = os.path.join(download_folder, "LIPIDMAPS")

# URLs
hmdb_url = "https://hmdb.ca/system/downloads/current/csf_metabolites.zip"
lm_url = "http://lipidmaps.org/files/?file=LMSD&ext=sdf.zip"

# Paths for downloaded ZIP files
hmdb_zip = os.path.join(hmdb_folder, "hmdb_metabolites.zip")
lm_zip = os.path.join(lm_folder, "LMSD.sdf.zip")

# ======= Version check config =======
UA = "checker/0.9 (+https://example.local)"
TIMEOUT = 30
TARGETS = {
    "HMDB": {
        "url": "https://hmdb.ca/downloads",
        "pattern": r"All Metabolites</td><td>(\d{4}-\d{2}-\d{2})</td><td><a data-toggle=\"modal\" data-target=\"#downloadModal\" data-whatever=\"/system/downloads/current/hmdb_metabolites.zip\"",
    },
    "LIPID MAPS": {
        "url": "https://lipidmaps.org/databases/lmsd/download",
        "pattern": r">LMSD\s+([\d-]+)\s+\(ZIP\)",
    },
}


def fetch_version(name: str, max_retries: int = 3) -> str:
    """Fetch remote version string using regex pattern, with retries."""
    target = TARGETS[name]

    for attempt in range(1, max_retries + 1):
        print(f"[{name}] Trying to fetch remote database version (attempt {attempt}/{max_retries})...")
        try:
            req = urllib.request.Request(target["url"], headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                charset = resp.headers.get_content_charset() or "utf-8"
                html = resp.read().decode(charset, errors="ignore")

            m = re.search(target["pattern"], html, flags=re.IGNORECASE)
            if not m:
                raise RuntimeError("Could not extract version from page")
            version = m.group(1)
            print(f"[{name}] Remote version found: {version}")
            return version

        except Exception as e:
            print(f"[{name}] Error fetching version: {e}")
            if attempt < max_retries:
                print(f"[{name}] Retrying...")
            else:
                print(f"[{name}] Failed after {max_retries} attempts.")
                raise


def download_zip(url, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    try:
        print(f"Downloading {os.path.basename(output_path)}...")
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()
        total_size = int(response.headers.get("Content-Length", 0))
        with (
            open(output_path, "wb") as f,
            tqdm(
                desc=f"Downloading {os.path.basename(output_path)}",
                total=total_size,
                unit="B",
                unit_scale=True,
            ) as bar,
        ):
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
                    bar.update(len(chunk))
        print(f"Downloaded {output_path}")
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {output_path}: {e}")
        retry_url = input("Provide a valid URL to retry: ").strip()
        if retry_url:
            download_zip(retry_url, output_path)


def extract_and_rename(zip_path, extract_to, base_filename, version):
    """Extract ZIP file, then rename main file to include version string."""
    try:
        print(f"Extracting {zip_path} to {extract_to}...")
        os.makedirs(extract_to, exist_ok=True)
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(extract_to)
            extracted_files = zip_ref.namelist()

        print(f"Extracted files: {extracted_files}")

        # Find first extracted file
        if extracted_files:
            extracted_file_path = os.path.join(extract_to, extracted_files[0])
            if not os.path.isfile(extracted_file_path):
                # Handle nested directory case
                for root, _, files in os.walk(extract_to):
                    for file in files:
                        extracted_file_path = os.path.join(root, file)
                        break

            target_file = f"{version}_{base_filename}"
            target_file_path = os.path.join(extract_to, target_file)
            shutil.move(extracted_file_path, target_file_path)
            print(f"Renamed to: {target_file_path}")
            return target_file_path
        else:
            print("No files extracted.")
            return None
    except Exception as e:
        print(f"Error extracting ZIP file: {e}")
        return None


if __name__ == "__main__":
    print("==== Getting current database versions ====")

    # Fetch remote versions
    try:
        hmdb_version = fetch_version("HMDB")
        lm_version = fetch_version("LIPID MAPS")
        print(f"HMDB remote version: {hmdb_version}")
        print(f"LIPID MAPS remote version: {lm_version}")
    except Exception as e:
        print(f"[ERROR] Could not fetch remote versions: {e}")
        sys.exit(1)

    print("==== Starting dataset download and extraction ====")
    # HMDB
    download_zip(hmdb_url, hmdb_zip)
    extract_and_rename(hmdb_zip, hmdb_folder, "hmdb_metabolites.xml", hmdb_version)

    # LIPID MAPS
    download_zip(lm_url, lm_zip)
    extract_and_rename(lm_zip, lm_folder, "structures.sdf", lm_version)

    print("==== Dataset download and extraction completed ====")
