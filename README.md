# HMDB Metabolite Explorer

A web-based interactive application built with **Shiny for Python** that enables exploration, filtering, and visualization of metabolite data.  
The tool focuses on **natural product classifications, lipid taxonomies, and biological properties**, providing an intuitive interface for researchers and analysts.  

The app now supports **automatic dataset updates** through the Admin Panel. Local dataset versions are extracted dynamically from the Parquet filename, and update checks compare them with the latest HMDB and LIPID MAPS releases. A dedicated test dataset is included for verifying the update workflow.


## License
This project is licensed under the terms of the [GNU General Public License v3.0](LICENSE).


## Features

- **Multi-criteria filtering**  
  Filter by pathway, tissue, lipid category, and more.  

- **Interactive visualizations**  
  Bar charts, histograms, and taxonomic breakdowns.  

- **Searchable DataGrid**  
  Styled and exportable table with rich filtering options.  

- **Detailed information panel**  
  Displays taxonomy and external database links for selected metabolites.  

- **Admin Panel**  
  Upload and preview new Parquet datasets  
  Run the enrichment pipeline for dataset updates  
  Check online database versions against local files  
  Dynamic loading of datasets with automatic version parsing  


## Recent Updates (2025-09-21)

- Added **update card** in Admin Panel to check online dataset versions  
- Integrated **run_pipeline.py** with Admin Panel for automated dataset updates  
- Enabled **dynamic parquet file loading** with version extraction from filenames  
- Added **test dataset** to validate Admin Panel update function  


## Usage

1. Prepare or update the required **Parquet file** using the enrichment pipeline.  
2. Run the application:  

   ```bash
   shiny run --reload app.py
   ```

3. Use the filter panel to refine results.  
4. Explore statistics, browse the DataGrid, and inspect detailed entries.  
5. Use the **Admin Panel** to upload datasets, run updates, or check for new HMDB and LIPID MAPS versions.  


## Requirements

Install the required dependencies before running the application:

```bash
pip install shiny shinywidgets pandas numpy plotly faicons pyarrow fastparquet
```

## Docker Usage

Clone the repository
```bash
git clone https://github.com/peterk140360/metabo-explorer.git
cd metabo-explorer
```

Build the image
```bash
docker build -t metabo-explorer:latest .
```

Run container
```bash
docker run -d --name metabo-explorer \
  -p 3838:3838 \
  -v /home/<user>/metabo-explorer/data:/app/data \
  --restart unless-stopped \
  metabo-explorer:latest
```

Check log
```bash
docker logs -f metabo-explorer
```

## To Do

- Freeze/fix first column in the DataGrid  
- Adopt the url in the collect-raw-data.py script from csf_ to hmdb_
- Implement console view in admin panel to show output of the pipeline run
- Check the remote server version error for updating the dataset


## Resources

- [Shiny for Python](https://shiny.posit.co/py/)  
- [Plotly](https://plotly.com/python/)  
- [Pandas](https://pandas.pydata.org/)  
- [HMDB](https://hmdb.ca/)  
- [LipidMaps](https://www.lipidmaps.org/)  
- [ChEBI](https://www.ebi.ac.uk/chebi/)  
- [NPClassifier](https://npclassifier.gnps2.org/)  
- [ClassyFire](http://classyfire.wishartlab.com/)  
