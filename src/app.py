###############################################################################
# Author:      Sebastian Peterka
# Company:     -
# File:        app.py
# Created:     2025-07-01
# Edited:      2025-07-31
#                           add 'not classified bar to barchart'
#                           change piechart to barchart
#                           delete "TOP {number}"" in barchart
#                           change order of plots
#                           add LM | Main Class as column in datagrid
#                           change name for taxonomy columns in datagrid
#                           no custom CSS for admin panel
#              2025-09-21  
#                           add update card in admin panel
#                           run pipeline for dataset update
#                           dynamic parquet file loading
#
# Description: The "HMDB Metabolite Explorer" is a web-based interactive
#              application built with Shiny for Python. It allows users to
#              explore, filter, and visualize metabolite data with a focus on
#              natural product classifications, lipid taxonomies, and biological
#              properties.
#
#              Key features include:
#                   * Multi-criteria filtering (pathway, tissue, lipid category, etc.)
#                   * Interactive visualizations (bar chart, histogram, pie chart)
#                   * Searchable and styled DataGrid with export functionality
#                   * Detailed information panel with taxonomy and external links
#                   * Admin panel with dataset management tools:
#                         - Update card to check for new dataset versions online
#                         - Trigger pipeline run for refreshing datasets
#                         - Dynamic parquet file loading and version parsing
#
# Usage:       1. Prepare or update the required Parquet file using the provided
#                 enrichment pipeline.
#              2. Run the app using:
#                     shiny run --reload app.py
#              3. Use the filter panel to narrow down results.
#              4. View statistics, browse the data grid, and inspect detailed entries.
#              5. Use the Admin Panel to upload new datasets, check update status,
#                 or re-run the pipeline.
#
# Comments:    * Requires the following Python packages:
#                   shiny, shinywidgets, pandas, numpy, plotly, faicons
#              * Ensure that the Parquet file naming scheme includes dataset
#                versions (e.g. 2021-11-17_hmdb_metabolites_classy_np_lm_2025-09-21_v5.parquet).
#              * The app extracts local dataset versions dynamically from the
#                filename for consistency with online checks.
#
# Todo:        * freeze/fix first column in the datagrid
#              * adopt the url in the collect-raw-data.py script from csf_ to hmdb_
#              * implement console view in admin panel to show output of the pipeline run
#
# Resources:   https://shiny.posit.co/py/
#              https://plotly.com/python/
#              https://pandas.pydata.org/
#              https://hmdb.ca/
#              https://www.lipidmaps.org/
#              https://www.ebi.ac.uk/chebi/
#              https://npclassifier.gnps2.org/
#              http://classyfire.wishartlab.com/
#
###############################################################################

# ---------------------------------------
# Imports
# ---------------------------------------
# Standard library
import re
import os
import sys
import pandas as pd
import numpy as np
import plotly.express as px
import subprocess
from faicons import icon_svg
from pathlib import Path
import urllib.request

# Shiny framework
from shiny import reactive
from shiny.express import input, render, ui
from shinywidgets import render_plotly


# ---------------------------------------
# Constants and Configuration
# ---------------------------------------
# Parquet file path dynamically
data_dir = Path(__file__).parent.parent / "data"
parquet_files = list(data_dir.glob("*.parquet"))

if not parquet_files:
    raise FileNotFoundError(f"No .parquet files found in {data_dir}")

# Pick the most recently modified parquet file
DATA_PATH = max(parquet_files, key=lambda f: f.stat().st_mtime)

print(f"Using Parquet file: {DATA_PATH}")

# DATA_PATH = (
#     # Path(__file__).parent.parent / "data" / "metabolites_enr_edit_np_lipidmaps_test.parquet"
#     Path(__file__).parent.parent / "data" / "metabolites_hmdb_np_lipidmaps.parquet"

# )

# Renaming dictionaries for display purposes
COLUMN_RENAMES_GRID = {
    # Top-level fields
    "accession": "HMDB ID",
    "status": "Status",
    "name": "Name",
    "chemical_formula": "Chemical Formula",
    "average_molecular_weight": "Avg. Molecular Weight",
    "smiles": "SMILES",
    "inchikey": "InChI Key",
    "pubchem_compound_id": "PubChem ID",
    "chebi_id": "ChEBI ID",
    "lm_id": "LipidMaps ID",
    # Flattened taxonomy fields
    "taxonomy_super_class": "HMDB | Superclass",
    # Flattened np_taxonomy fields
    "np_taxonomy_pathway": "NP Classifier | Pathway",
    # Flattened lm_taxonomy fields
    "lm_taxonomy_category": "Lipid Maps | Category",
    "lm_taxonomy_main_class": "Lipid Maps | Main Class",
    # Flattened biological_properties fields
    "biological_properties_cellular_locations": "Cellular Locations",
    "biological_properties_biospecimen_locations": "Biospecimen Locations",
    "biological_properties_tissue_locations": "Tissue Locations",
}
COLUMN_RENAMES_EXPORT = {
    # Top-level fields
    "accession": "HMDB_ID",
    "status": "Status",
    "name": "Name",
    "description": "Description",
    "chemical_formula": "Chemical_Formula",
    "average_molecular_weight": "Avg_Molecular_Weight",
    "iupac_name": "IUPAC_Name",
    "smiles": "SMILES",
    "inchikey": "InChIKey",
    "pubchem_compound_id": "PubChem_ID",
    "chebi_id": "ChEBI_ID",
    "wikipedia_id": "Wikipedia_ID",
    "lm_id": "LipidMaps_ID",
    # Flattened taxonomy fields
    "taxonomy_kingdom": "HMDB_Kingdom",
    "taxonomy_super_class": "HMDB_Superclass ",
    "taxonomy_class": "HMDB_Class",
    "taxonomy_sub_class": "HMDB_Subclass",
    "taxonomy_direct_parent": "HMDB_Direct_Parent",
    # Flattened np_taxonomy fields
    "np_taxonomy_pathway": "NP_Pathway",
    "np_taxonomy_super_class": "NP_Superclass",
    "np_taxonomy_class": "NP_Class",
    # Flattened lm_taxonomy fields
    "lm_taxonomy_category": "LM_Category",
    "lm_taxonomy_main_class": "LM_Main_Class",
    "lm_taxonomy_sub_class": "LM_Sub_Class",
    # Flattened biological_properties fields
    "biological_properties_cellular_locations": "Cellular_Locations",
    "biological_properties_biospecimen_locations": "Biospecimen_Locations",
    "biological_properties_tissue_locations": "Tissue_Locations",
}
COLUMN_RENAMES_DETAIL = {
    # Top-level fields
    "accession": "HMDB ID",
    "status": "Status",
    "name": "Name",
    "description": "Description",
    "chemical_formula": "Chemical Formula",
    "average_molecular_weight": "Avg. Molecular Weight",
    "iupac_name": "IUPAC Name",
    "smiles": "SMILES",
    "inchikey": "InChI Key",
    "pubchem_compound_id": "PubChem ID",
    "chebi_id": "ChEBI ID",
    "wikipedia_id": "Wikipedia ID",
    "lm_id": "LipidMaps ID",
    # Flattened taxonomy fields
    "taxonomy_kingdom": "Kingdom",
    "taxonomy_super_class": "Superclass ",
    "taxonomy_class": "Class ",
    "taxonomy_sub_class": "Subclass",
    "taxonomy_direct_parent": "Direct Parent",
    # Flattened np_taxonomy fields
    "np_taxonomy_pathway": "Pathway",
    "np_taxonomy_super_class": "Superclass",
    "np_taxonomy_class": "Class",
    # Flattened lm_taxonomy fields
    "lm_taxonomy_category": "Category",
    "lm_taxonomy_main_class": "Main-Class",
    "lm_taxonomy_sub_class": "Sub-Class",
    # Flattened biological_properties fields
    "biological_properties_cellular_locations": "Cellular Locations",
    "biological_properties_biospecimen_locations": "Biospecimen Locations",
    "biological_properties_tissue_locations": "Tissue Locations",
}

# Columns to hide in UI
HIDDEN_COLUMNS = [
    "status",
    "description",
    "iupac_name",
    "wikipedia_id",
    "taxonomy_kingdom",
    "taxonomy_class",
    "taxonomy_sub_class",
    "taxonomy_direct_parent",
    "np_taxonomy_super_class",
    "np_taxonomy_class",
    "lm_taxonomy_sub_class",
]

# Columns that represent biological properties and require list normalization
BIO_PROP = [
    "biological_properties_cellular_locations",
    "biological_properties_biospecimen_locations",
    "biological_properties_tissue_locations",
]


# ---------------------------------------
# Main Application Logic
# ---------------------------------------
# Load metabolite data from Parquet file
df_default = pd.read_parquet(DATA_PATH)


# ---------------------------------------
# Reactive Expressions
# ---------------------------------------
@reactive.calc
def dataset() -> pd.DataFrame:
    """
    Return the current dataset, using the uploaded file if selected.
    Falls back to the default .parquet file if no upload is used.
    """

    def normalize_biological_columns(df: pd.DataFrame) -> pd.DataFrame:
        """
        # Convert biological property columns to plain lists
        # Example: numpy.ndarray ['Blood' 'Urine'] → list ['Blood', 'Urine']
        """
        for col in BIO_PROP:
            if col in df.columns:
                df[col] = df[col].apply(
                    lambda x: x.tolist() if isinstance(x, np.ndarray) else x
                )
        return df

    file_info = input.parquet_file()

    if input.upload_btn() == 0 or not file_info:
        return normalize_biological_columns(df_default.copy())

    try:
        df_uploaded = pd.read_parquet(file_info[0]["datapath"])
        return normalize_biological_columns(df_uploaded)
    except Exception:
        return normalize_biological_columns(df_default.copy())


@reactive.calc
def filtered_view() -> pd.DataFrame:
    """
    Return a filtered copy of the dataset based on current input selections.
    Each filter is applied independently if a corresponding input is set.
    """
    df_filter = dataset().copy()

    # Filter by NP pathway classification
    if input.pathway_filter():
        df_filter = df_filter[
            df_filter["np_taxonomy_pathway"].isin(input.pathway_filter())
        ]

    # Filter by cellular locations (list column)
    if input.cell_location_filter():
        df_filter = df_filter[
            df_filter["biological_properties_cellular_locations"].apply(
                lambda x: any(loc in x for loc in input.cell_location_filter())
                if isinstance(x, list)
                else False
            )
        ]

    # Filter by lipid category classification
    if input.lipid_category_filter():
        df_filter = df_filter[
            df_filter["lm_taxonomy_category"].isin(input.lipid_category_filter())
        ]

    # Filter by biospecimen locations (list column)
    if input.biospecimen_location_filter():
        df_filter = df_filter[
            df_filter["biological_properties_biospecimen_locations"].apply(
                lambda x: any(loc in x for loc in input.biospecimen_location_filter())
                if isinstance(x, list)
                else False
            )
        ]

    # Filter by superclass taxonomy
    if input.superclass_filter():
        df_filter = df_filter[
            df_filter["taxonomy_super_class"].isin(input.superclass_filter())
        ]

    # Filter by tissue locations (list column)
    if input.tissue_location_filter():
        df_filter = df_filter[
            df_filter["biological_properties_tissue_locations"].apply(
                lambda x: any(loc in x for loc in input.tissue_location_filter())
                if isinstance(x, list)
                else False
            )
        ]

    return df_filter


# ---------------------------------------
# UI Setup
# ---------------------------------------
# Apply custom CSS for scrollable selectize UI
ui.tags.style("""
.selectize-scroll .selectize-control.multi .selectize-input {
    max-height: 65px;
    overflow-y: auto !important;
    overflow-x: hidden;
    flex-wrap: nowrap !important;
}

/* Force table headers in admin panel to left align */
th { text-align: left !important; }
td { text-align: left; }
""")


# ---------------------------------------
# Define UI Header
# ---------------------------------------
ui.div(
    ui.h4("HMDB Metabolite Explorer", class_="mb-0"),
    class_="bg-primary text-white p-3 mb-3 rounded",
    style="position: static; top: 0; z-index: 1000;",
)


# ---------------------------------------
# Tabs: Main View + Admin Panel
# ---------------------------------------
with ui.navset_pill(id="main_tab"):
    # --------- Main View Tab ----------
    with ui.nav_panel("Main View", value="main"):
        # ---------------------------------------
        # UI TOP ROW: with Filters and Plot
        # ---------------------------------------
        with ui.layout_columns(col_widths=(6, 6)):
            # Main filter card
            with ui.card():
                ui.card_header("Filter")

                with ui.layout_column_wrap(width=1 / 2):
                    # Pathway filter (NP-Classifier)
                    with ui.div(style="height: 130px; margin-bottom: 10px;"):
                        with ui.panel_well(class_="h-100 selectize-scroll"):
                            ui.p("Pathway Filter (NP-Classifier)")
                            ui.input_selectize(
                                id="pathway_filter",
                                label="",
                                choices={
                                    "NP-Classifier Pathways": {
                                        "Fatty acids": "Fatty acids",
                                        "Polyketides": "Polyketides",
                                        "Amino acids and Peptides": "Amino acids and Peptides",
                                        "Terpenoids": "Terpenoids",
                                        "Shikimates and Phenylpropanoids": "Shikimates and Phenylpropanoids",
                                        "Alkaloids": "Alkaloids",
                                        "Carbohydrates": "Carbohydrates",
                                    }
                                },
                                multiple=True,
                                options={
                                    "placeholder": "Select pathways...",
                                    "plugins": ["clear_button"],
                                    "render": ui.js_eval(
                                        '{"optgroup_header": function(data, e) {return "<div style=\'font-weight:bold;\'>" + e(data.label) + "</div>";}, "option": function(item, e) {return "<div style=\'padding-left: 1rem;\'>" + e(item.label) + "</div>";}}'
                                    ),
                                },
                            )

                    # Cellular location filter
                    with ui.div(style="height: 130px; margin-bottom: 10px;"):
                        with ui.panel_well(class_="h-100 selectize-scroll"):
                            ui.p("Cellular Location Filter")
                            ui.input_selectize(
                                id="cell_location_filter",
                                label="",
                                choices={
                                    "HMDB Cellular Locations": {
                                        "Extracellular": "Extracellular",
                                        "Lysosome": "Lysosome",
                                        "Mitochondria": "Mitochondria",
                                        "Nucleus": "Nucleus",
                                        "Membrane": "Membrane",
                                        "Cytoplasm": "Cytoplasm",
                                        "Endoplasmic reticulum": "Endoplasmic reticulum",
                                    }
                                },
                                multiple=True,
                                options={
                                    "placeholder": "Select cellular locations...",
                                    "plugins": ["clear_button"],
                                    "render": ui.js_eval(
                                        '{"optgroup_header": function(data, e) {return "<div style=\'font-weight:bold;\'>" + e(data.label) + "</div>";}, "option": function(item, e) {return "<div style=\'padding-left: 1rem;\'>" + e(item.label) + "</div>";}}'
                                    ),
                                },
                            )

                    # Lipid category filter (Lipid Maps)
                    with ui.div(style="height: 130px; margin-bottom: 10px;"):
                        with ui.panel_well(class_="h-100 selectize-scroll"):
                            ui.p("Category Filter (LipidMaps)")
                            ui.input_selectize(
                                id="lipid_category_filter",
                                label="",
                                choices={
                                    "Lipid Maps Categories": {
                                        "Fatty Acyls [FA]": "Fatty Acyls",
                                        "Glycerolipids [GL]": "Glycerolipids",
                                        "Glycerophospholipids [GP]": "Glycerophospholipids",
                                        "Polyketides [PK]": "Polyketides",
                                        "Prenol Lipids [PR]": "Prenol Lipids",
                                        "Saccharolipids [SL]": "Saccharolipids",
                                        "Sphingolipids [SP]": "Sphingolipids",
                                        "Sterol Lipids [ST]": "Sterol Lipids",
                                    }
                                },
                                multiple=True,
                                options={
                                    "placeholder": "Select lipid categories...",
                                    "plugins": ["clear_button"],
                                    "render": ui.js_eval(
                                        '{"optgroup_header": function(data, e) {return "<div style=\'font-weight:bold;\'>" + e(data.label) + "</div>";}, "option": function(item, e) {return "<div style=\'padding-left: 1rem;\'>" + e(item.label) + "</div>";}}'
                                    ),
                                },
                            )

                    # Biospecimen location filter
                    with ui.div(style="height: 130px; margin-bottom: 10px;"):
                        with ui.panel_well(class_="h-100 selectize-scroll"):
                            ui.p("Biospecimen Location Filter")
                            ui.input_selectize(
                                id="biospecimen_location_filter",
                                label="",
                                choices={
                                    "HMDB Biospecimen Locations": {
                                        "Blood": "Blood",
                                        "Cerebrospinal Fluid (CSF)": "Cerebrospinal Fluid (CSF)",
                                        "Feces": "Feces",
                                        "Saliva": "Saliva",
                                        "Urine": "Urine",
                                    }
                                },
                                multiple=True,
                                options={
                                    "placeholder": "Select biospecimen locations...",
                                    "plugins": ["clear_button"],
                                    "render": ui.js_eval(
                                        '{"optgroup_header": function(data, e) {return "<div style=\'font-weight:bold;\'>" + e(data.label) + "</div>";}, "option": function(item, e) {return "<div style=\'padding-left: 1rem;\'>" + e(item.label) + "</div>";}}'
                                    ),
                                },
                            )

                    # Superclass filter (ClassyFire)
                    with ui.div(style="height: 130px; margin-bottom: 10px;"):
                        with ui.panel_well(class_="h-100 selectize-scroll"):
                            ui.p("Superclass Filter (ClassyFire)")
                            ui.input_selectize(
                                id="superclass_filter",
                                label="",
                                choices={
                                    "HMDB Superclass": {
                                        "Acetylides": "Acetylides",
                                        "Alkaloids and derivatives": "Alkaloids and derivatives",
                                        "Allenes": "Allenes",
                                        "Benzenoids": "Benzenoids",
                                        "Carbenes": "Carbenes",
                                        "Hydrocarbon derivatives": "Hydrocarbon derivatives",
                                        "Hydrocarbons": "Hydrocarbons",
                                        "Lignans, neolignans and related compounds": "Lignans, neolignans and related compounds",
                                        "Lipids and lipid-like molecules": "Lipids and lipid-like molecules",
                                        "Nucleosides, nucleotides, and analogues": "Nucleosides, nucleotides, and analogues",
                                        "Organic 1,3-dipolar compounds": "Organic 1,3-dipolar compounds",
                                        "Organic Polymers": "Organic Polymers",
                                        "Organic acids and derivatives": "Organic acids and derivatives",
                                        "Organic anions": "Organic anions",
                                        "Organic cations": "Organic cations",
                                        "Organic nitrogen compounds": "Organic nitrogen compounds",
                                        "Organic oxygen compounds": "Organic oxygen compounds",
                                        "Organic salts": "Organic salts",
                                        "Organic zwitterions": "Organic zwitterions",
                                        "Organohalogen compounds": "Organohalogen compounds",
                                        "Organoheterocyclic compounds": "Organoheterocyclic compounds",
                                        "Organometallic compounds": "Organometallic compounds",
                                        "Organophosphorus compounds": "Organophosphorus compounds",
                                        "Organopnictogen compounds": "Organopnictogen compounds",
                                        "Organosulfur compounds": "Organosulfur compounds",
                                        "Phenylpropanoids and polyketides": "Phenylpropanoids and polyketides",
                                    }
                                },
                                multiple=True,
                                options={
                                    "placeholder": "Select superclasses...",
                                    "plugins": ["clear_button"],
                                    "render": ui.js_eval(
                                        '{"optgroup_header": function(data, e) {return "<div style=\'font-weight:bold;\'>" + e(data.label) + "</div>";}, "option": function(item, e) {return "<div style=\'padding-left: 1rem;\'>" + e(item.label) + "</div>";}}'
                                    ),
                                },
                            )

                    # Tissue location filter
                    with ui.div(style="height: 130px; margin-bottom: 10px;"):
                        with ui.panel_well(class_="h-100 selectize-scroll"):
                            ui.p("Tissue Location Filter")
                            ui.input_selectize(
                                id="tissue_location_filter",
                                label="",
                                choices={
                                    "HMDB Tissue Locations": {
                                        "Adrenal Gland": "Adrenal Gland",
                                        "Brain": "Brain",
                                        "Epidermis": "Epidermis",
                                        "Fibroblasts": "Fibroblasts",
                                        "Liver": "Liver",
                                        "Testis": "Testis",
                                    }
                                },
                                multiple=True,
                                options={
                                    "placeholder": "Select tissue locations...",
                                    "plugins": ["clear_button"],
                                    "render": ui.js_eval(
                                        '{"optgroup_header": function(data, e) {return "<div style=\'font-weight:bold;\'>" + e(data.label) + "</div>";}, "option": function(item, e) {return "<div style=\'padding-left: 1rem;\'>" + e(item.label) + "</div>";}}'
                                    ),
                                },
                            )

            # Plot card
            with ui.card(full_screen=True):
                with ui.card_header(
                    class_="d-flex justify-content-between align-items-center"
                ):
                    "Summary Statistics & Plots"

                    # Popover with plot type selector
                    with ui.popover(title="Select plot type", placement="top"):
                        icon_svg("ellipsis")
                        ui.input_select(
                            "plot_type",
                            label=None,
                            choices={
                                "pathway_bar": "NP Pathways",
                                "lipid_bar": "Lipid Categorys",
                                "mol_weight_hist": "Molecular Weight Distribution",
                            },
                        )

                    def value_counts_with_unclassified(series: pd.Series) -> pd.Series:
                        """
                        Count values in a Series and include an 'Not Classified' bar for empty or NaN entries.

                        Args:
                            series: The input Series to count values from.

                        Returns:
                            A Series of value counts
                        """
                        s_clean = series.fillna("").astype(str).str.strip()
                        counts = s_clean[s_clean != ""].value_counts()
                        n_unclassified = (s_clean == "").sum()
                        if n_unclassified > 0:
                            counts["Not Classified"] = n_unclassified
                        return counts

                @render_plotly
                def summary_plot():
                    df_plot = filtered_view()

                    # Return placeholder if no data matches the filters
                    if df_plot.empty:
                        fig = px.scatter()
                        fig.update_layout(
                            xaxis=dict(visible=False),
                            yaxis=dict(visible=False),
                            annotations=[
                                dict(
                                    text="No data to show",
                                    x=0.5,
                                    y=0.5,
                                    xref="paper",
                                    yref="paper",
                                    showarrow=False,
                                    font=dict(size=18),
                                )
                            ],
                        )
                        return fig

                    plot_type = input.plot_type()
                    color = "#007bc2"

                    # Bar chart: NP pathways
                    if plot_type == "pathway_bar":
                        counts = value_counts_with_unclassified(
                            df_plot["np_taxonomy_pathway"]
                        )
                        n = counts.sum()

                        fig = px.bar(
                            x=counts.index,
                            y=counts.values,
                            labels={"x": "Pathway", "y": "Count"},
                            title=f"NP Classifier Pathways (n = {n})",
                            color_discrete_sequence=[color],
                        )
                        fig.update_layout(xaxis_tickangle=-45)

                    # Bar chart: Lipid categories
                    elif plot_type == "lipid_bar":
                        counts = value_counts_with_unclassified(
                            df_plot["lm_taxonomy_category"]
                        )
                        n = counts.sum()

                        fig = px.bar(
                            x=counts.index,
                            y=counts.values,
                            labels={"x": "Lipid Category", "y": "Count"},
                            title=f"Lipid Categories (n = {n})",
                            color_discrete_sequence=[color],
                        )
                        fig.update_layout(xaxis_tickangle=-45)

                    # Histogram: Molecular weight distribution
                    elif plot_type == "mol_weight_hist":
                        df_plot["average_molecular_weight"] = pd.to_numeric(
                            df_plot["average_molecular_weight"], errors="coerce"
                        )

                        # Count non-empty, non-NaN values in the average_molecular_weight column
                        n = df_plot["average_molecular_weight"].count()

                        fig = px.histogram(
                            df_plot.dropna(subset=["average_molecular_weight"]),
                            x="average_molecular_weight",
                            nbins=30,
                            title=f"Molecular Weight Distribution (n = {n})",
                            color_discrete_sequence=[color],
                        )
                        fig.update_layout(
                            xaxis_title="Avg. Molecular Weight", yaxis_title="Count"
                        )

                    return fig

        # ---------------------------------------
        # UI MIDDLE ROW: DataGrid
        # ---------------------------------------
        with ui.card(style="height: 600px; overflow-y: auto;", full_screen=True):

            @render.data_frame
            def metabolites_table():
                df_datagrid = filtered_view().copy()

                # Convert list-type columns to comma-separated strings
                for col in df_datagrid.columns:
                    if df_datagrid[col].apply(lambda x: isinstance(x, list)).any():
                        df_datagrid[col] = df_datagrid[col].apply(
                            lambda x: ", ".join(x) if isinstance(x, list) else x
                        )

                # Remove columns hidden from the UI
                df_datagrid.drop(columns=HIDDEN_COLUMNS, inplace=True, errors="ignore")

                # Rename columns for display
                df_datagrid.rename(columns=COLUMN_RENAMES_GRID, inplace=True)

                # Define column-specific styling for the datagrid
                styles = []
                for i, col in enumerate(df_datagrid.columns):
                    style = {
                        "cols": [i],
                        "style": {
                            "minWidth": "50px",
                            "maxWidth": "300px",
                            "overflow": "hidden",
                            "text-overflow": "ellipsis",
                            "white-space": "nowrap",
                        },
                    }

                    # Align selected columns
                    if col == "Avg. Molecular Weight":
                        style["class"] = "text-center"
                    elif col in ["PubChem ID", "ChEBI ID"]:
                        style["class"] = "text-end"

                    styles.append(style)

                return render.DataGrid(
                    df_datagrid,
                    filters=True,
                    styles=styles,
                    summary="Showing rows {start} to {end} of {total}",
                    selection_mode="row",
                )

            @render.download(label="Download CSV")
            def download():
                """
                This download handler dynamically creates a CSV file from the filtered dataset
                and returns its path for download.
                """

                # Create a copy of the current filtered dataset
                df_export = filtered_view().copy()

                # Define the path where the CSV will be saved
                export_path = Path(__file__).parent / "filtered_export.csv"

                # Convert list-type columns to comma-separated strings for cleaner export
                for col in df_export.columns:
                    if df_export[col].apply(lambda x: isinstance(x, list)).any():
                        df_export[col] = df_export[col].apply(
                            lambda x: ", ".join(x) if isinstance(x, list) else x
                        )

                # Rename columns according to predefined export mapping
                df_export.rename(columns=COLUMN_RENAMES_EXPORT, inplace=True)

                # Write the DataFrame to a CSV file without the index
                df_export.to_csv(export_path, index=False)

                # Return the file path so Shiny Express can serve it as a download
                return str(export_path)

        # ---------------------------------------
        # UI BOTTOM ROW: Detail Panel
        # ---------------------------------------
        with ui.card(style="height: 500px; overflow-y: auto;", full_screen=True):

            @render.ui
            def selected_row_details():
                df_display = filtered_view().copy()

                # Rename columns for detailed display
                df_display.rename(columns=COLUMN_RENAMES_DETAIL, inplace=True)

                # Get the selected row from the DataGrid
                selection = metabolites_table.cell_selection()
                selected_rows = selection["rows"]

                # Show placeholder if no row is selected or the dataset is empty
                if not selected_rows or df_display.empty:
                    return ui.div(
                        "Select a row to see full details.", class_="text-muted mt-3"
                    )

                # Extract the selected row
                index = df_display.index[selected_rows[0]]
                row = df_display.loc[index]
                table_rows = []

                # ---------------------------
                # Helper: Determine indentation level for nested/hierarchical fields
                # ---------------------------
                def get_indent_level(label: str) -> int:
                    for prefix, level in [
                        ("Kingdom", 1),
                        ("Superclass ", 2),
                        ("Class ", 3),
                        ("Subclass", 4),
                        ("Direct Parent", 5),
                        ("Pathway", 1),
                        ("Superclass", 2),
                        ("Class", 3),
                        ("Category", 1),
                        ("Main-Class", 2),
                        ("Sub-Class", 3),
                        ("Cellular Locations", 1),
                        ("Biospecimen Locations", 1),
                        ("Tissue Locations", 1),
                    ]:
                        if label.startswith(prefix):
                            return level
                    return 0

                # ---------------------------
                # Helper: Render section header row
                # ---------------------------
                def section_header(title: str):
                    return ui.tags.tr(
                        ui.tags.td(
                            ui.tags.strong(title),
                            colspan=2,
                            style="background-color: #f0f8ff;",
                        )
                    )

                # Flags to avoid repeating section headers
                inserted_sections = {
                    "HMDB": False,
                    "NP": False,
                    "LipidMaps": False,
                    "Biological": False,
                }

                # Render each field in the selected row
                for col, val in row.items():
                    # Format list values as comma-separated strings
                    if isinstance(val, list):
                        val = ", ".join(val)

                    # Render known ID fields as clickable external links
                    if col == "HMDB ID" and pd.notna(val):
                        val = ui.a(
                            val,
                            href=f"https://hmdb.ca/metabolites/{val}",
                            target="_blank",
                        )
                    elif col == "LipidMaps ID" and pd.notna(val):
                        val = ui.a(
                            val,
                            href=f"https://www.lipidmaps.org/databases/lmsd/{val}",
                            target="_blank",
                        )
                    elif col == "PubChem ID" and pd.notna(val):
                        val = ui.a(
                            val,
                            href=f"https://pubchem.ncbi.nlm.nih.gov/compound/{val}",
                            target="_blank",
                        )
                    elif col == "ChEBI ID" and pd.notna(val):
                        chebi = str(val).replace("CHEBI:", "").upper()
                        val = ui.a(
                            f"CHEBI:{chebi}",
                            href=f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{chebi}",
                            target="_blank",
                        )

                    # Insert taxonomy/section headers before hierarchical groups
                    if col.startswith("Kingdom") and not inserted_sections["HMDB"]:
                        table_rows.append(section_header("HMDB Taxonomy"))
                        inserted_sections["HMDB"] = True
                    elif col.startswith("Pathway") and not inserted_sections["NP"]:
                        table_rows.append(section_header("NP Taxonomy"))
                        inserted_sections["NP"] = True
                    elif (
                        col.startswith("Category")
                        and not inserted_sections["LipidMaps"]
                    ):
                        table_rows.append(section_header("Lipid Maps Taxonomy"))
                        inserted_sections["LipidMaps"] = True
                    elif (
                        col == "Cellular Locations"
                        and not inserted_sections["Biological"]
                    ):
                        table_rows.append(section_header("Biological Properties"))
                        inserted_sections["Biological"] = True

                    # Determine indentation for hierarchical fields
                    indent_px = get_indent_level(col) * 20
                    style_val = f"padding: 6px 12px; padding-left: {indent_px}px;"
                    if col == "Description":
                        style_val += " text-align: justify;"

                    # Final table row layout
                    table_rows.append(
                        ui.tags.tr(
                            ui.tags.td(
                                col,
                                style=(
                                    f"vertical-align: top; padding: 6px 12px; padding-left: {indent_px}px; width: 30%;"
                                    + (
                                        ""
                                        if any(
                                            col.startswith(prefix)
                                            for prefix in [
                                                "Kingdom",
                                                "Superclass ",
                                                "Class ",
                                                "Subclass",
                                                "Direct Parent",
                                                "Pathway",
                                                "Superclass",
                                                "Class",
                                                "Category",
                                                "Main-Class",
                                                "Sub-Class",
                                                "Cellular Locations",
                                                "Biospecimen Locations",
                                                "Tissue Locations",
                                            ]
                                        )
                                        else " font-weight: bold;"
                                    )
                                ),
                            ),
                            ui.tags.td(
                                ui.div(val, style="text-align: justify;")
                                if col == "Description"
                                else val,
                                style=style_val,
                            ),
                        )
                    )

                # Determine metabolite name for header (fallback to ID)
                metabolite_name = (
                    row.get("Name") or row.get("HMDB ID") or "Selected Metabolite"
                )

                # Final output: card title and table of details
                return ui.div(
                    ui.div(
                        ui.tags.h5(
                            f"Full Details for {metabolite_name}", class_="mb-0"
                        ),
                        class_="p-2",
                        style="background-color: #e6f2ff; border-radius: 5px;",
                    ),
                    ui.tags.table(
                        *table_rows,
                        style="width: 100%; border-collapse: collapse; font-size: 0.95rem;",
                    ),
                    class_="mt-2",
                )

    # --------- Admin Panel Tab (via menu) ----------
    with ui.nav_menu(title="", icon=icon_svg("ellipsis-vertical"), align="right"):
        with ui.nav_panel("Admin Panel", value="admin"):
            ui.h3("Admin Panel")
            ui.p("This panel is for administrators only.")

            with ui.card():
                ui.card_header([icon_svg("arrows-rotate"), " Check for Updates"])
                ui.input_action_button(
                    "check_updates", "Check for update", class_="btn-warning"
                )

                outdated_flag = reactive.Value(False)

                @render.table
                @reactive.event(input.check_updates)
                def update_status():
                    UA = "checker/0.8 (+https://example.local)"
                    TIMEOUT = 30

                    # Get local versions dynamically from the parquet filename
                    parquet_path = str(DATA_PATH)  # your global Path object
                    filename = os.path.basename(parquet_path)

                    m = re.match(
                        r"(\d{4}-\d{2}-\d{2})_hmdb_metabolites.*_lm_(\d{4}-\d{2}-\d{2})",
                        filename,
                    )
                    if m:
                        hmdb_version, lm_version = m.group(1), m.group(2)
                    else:
                        hmdb_version, lm_version = "unknown", "unknown"

                    LOCAL_VERSIONS = {
                        "HMDB": hmdb_version,
                        "LIPID MAPS": lm_version,
                    }

                    # ---- Regex patterns to scrape remote versions ----
                    TARGETS = [
                        {
                            "name": "HMDB",
                            "url": "https://hmdb.ca/downloads",
                            "pattern": r"All Metabolites</td><td>(\d{4}-\d{2}-\d{2})</td><td><a data-toggle=\"modal\" data-target=\"#downloadModal\" data-whatever=\"/system/downloads/current/hmdb_metabolites.zip\"",
                        },
                        {
                            "name": "LIPID MAPS",
                            "url": "https://lipidmaps.org/databases/lmsd/download",
                            "pattern": r">LMSD\s+([\d-]+)\s+\(ZIP\)",
                        },
                    ]

                    def fetch(url: str) -> str:
                        req = urllib.request.Request(url, headers={"User-Agent": UA})
                        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                            charset = resp.headers.get_content_charset() or "utf-8"
                            return resp.read().decode(charset, errors="ignore")

                    results = []
                    any_outdated = False
                    for t in TARGETS:
                        local_version = LOCAL_VERSIONS.get(t["name"], "n/a")
                        try:
                            html = fetch(t["url"])
                            m = re.search(t["pattern"], html, flags=re.IGNORECASE)
                            if m:
                                remote_version = m.group(1)
                                status = (
                                    "UP TO DATE"
                                    if local_version == remote_version
                                    else "OUTDATED"
                                )
                                if status == "OUTDATED":
                                    any_outdated = True
                                results.append(
                                    {
                                        "Database": t["name"],
                                        "Local Version": local_version,
                                        "Remote Version": remote_version,
                                        "Status": status,
                                    }
                                )
                            else:
                                results.append(
                                    {
                                        "Database": t["name"],
                                        "Local Version": local_version,
                                        "Remote Version": "not found",
                                        "Status": "ERROR: could not extract version",
                                    }
                                )
                                any_outdated = True
                        except Exception as e:
                            results.append(
                                {
                                    "Database": t["name"],
                                    "Local Version": local_version,
                                    "Remote Version": "error",
                                    "Status": f"ERROR: {e}",
                                }
                            )
                            any_outdated = True

                    outdated_flag.set(any_outdated)
                    return pd.DataFrame(results)

                # ---- Password + Update button (only if outdated) ----
                @render.ui
                def update_controls():
                    if outdated_flag.get():
                        return ui.div(
                            ui.input_password("admin_pw", "Enter password:"),
                            ui.input_action_button(
                                "update_db", "Update Database", class_="btn-danger mt-2"
                            ),
                        )
                    return None

                # ---- Run pipeline when button pressed ----
                @reactive.effect
                @reactive.event(input.update_db)
                def _():
                    pw = input.admin_pw()
                    if pw != "admin":
                        ui.modal_show(
                            ui.modal(
                                "Incorrect password!", title="Error", easy_close=True
                            )
                        )
                        return

                    script_path = Path(__file__).parent / "run_pipeline.py"
                    print("Launching pipeline at:", script_path)

                    try:
                        subprocess.Popen([sys.executable, str(script_path)])
                        ui.modal_show(
                            ui.modal(
                                "Pipeline started successfully!",
                                title="Success",
                                easy_close=True,
                            )
                        )
                    except Exception as e:
                        ui.modal_show(
                            ui.modal(
                                f"Failed to start pipeline: {e}",
                                title="Error",
                                easy_close=True,
                            )
                        )

            with ui.card():
                ui.card_header([icon_svg("github"), " GitHub Repository"])
                ui.p(
                    [
                        "To manually update dataset visit the GitHub repository: ",
                        ui.a(
                            "GitHub repository",
                            href="https://github.com/peterk140360/Metabo-explorer",
                            target="_blank",
                        ),
                        ".",
                    ]
                )

            with ui.card(style="height: 200px"):
                ui.card_header(
                    [
                        icon_svg("upload"),
                        "Upload a new .parquet file to inspect its structure",
                    ]
                )
                ui.input_file(
                    id="parquet_file",
                    label="Choose a .parquet file",
                    accept=[".parquet"],
                    multiple=False,
                    width="100%",
                    placeholder="No file selected",
                )

            with ui.card(style="height: 300px; overflow-y: auto;"):
                ui.card_header([icon_svg("eye"), " Preview Summary"])

                @render.table
                def parquet_summary():
                    file_info = input.parquet_file()
                    if file_info is None:
                        return pd.DataFrame({"Message": ["No file uploaded."]})

                    try:
                        df = pd.read_parquet(file_info[0]["datapath"])
                    except Exception as e:
                        return pd.DataFrame({"Error": [str(e)]})

                    row_count = df.shape[0]
                    column_count = df.shape[1]
                    column_names = ", ".join(df.columns)

                    return pd.DataFrame(
                        {
                            "Row Count": [row_count],
                            "Column Count": [column_count],
                            "Column Names": [column_names],
                        }
                    )

            with ui.card():
                ui.input_action_button(
                    "upload_btn",
                    "Use uploaded file",
                    icon=icon_svg("upload"),
                    class_="btn-primary",
                )

                @render.ui
                def upload_status():
                    return ui.tags.div(
                        "File uploaded successfully!", class_="text-success mt-2"
                    )


# Define UI Footer
# ---------------------------------------
ui.div(
    ui.p(
        "© 2025 HMDB Metabolite Explorer — Built with Shiny for Python", class_="mb-0"
    ),
    class_="bg-light text-center text-muted p-3 mt-4 rounded",
    style="position: static; bottom: 0; z-index: 1000;",
)
