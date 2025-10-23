# TC-Hazard-Modelling-Leeward-Islands

This repository contains the Python scripts and Jupyter notebooks used for modelling Tropical Cyclone (TC) Wind, Storm Surge, and Rainfall hazards for a set of islands in the Leeward Islands region, primarily using the CLIMADA hazard modelling platform. The workflow is concentrated in a set of functions in the TC_hazard_modelling.py script and is further demonstrated across five Jupyter notebooks.

### Repository Structure Overview
The repository is organised to separate input data, pre-processing files, output results, and the core code logic.

| Name                        | Type              | Description                                                                 |
|-----------------------------|-------------------|-----------------------------------------------------------------------------|
| 0-level-Leeward++           | Folder            | Contains the initial raw data or master file (e.g., track data).             |
| ERA5_data                   | Folder            | Output folder for raw ERA5 NetCDF files downloaded by the CDS API.           |
| ERA5_parameters             | Folder            | Output folder for extracted ERA5 parameters (T600, U850, V850), saved as CSV.|
| rain_output                 | Folder            | Output folder for the final rainfall hazard results (.csv).                  |
| SRTM15Plus                  | Folder            | Input folder for DEM data (GeoTIFF files) used in storm surge modeling.      |
| STORM_PresentDayClimate     | Folder            | Input folder for synthetic storm tracks (present-day climate).                |
| surge_figures               | Folder            | Output folder for generated storm surge plots.                               |
| surge_output                | Folder            | Output folder for the final storm surge hazard results (.csv).                |
| wind_figures                | Folder            | Output folder for generated wind hazard plots.                               |
| wind_output                 | Folder            | Output folder for the final wind hazard results (.csv).                      |
| TC_hazard_modelling.py      | Python Script     | Core code file containing all hazard modeling functions.                     |
| IBTrACS_precipitation.ipynb | Jupyter Notebook  | Example workflow for rainfall hazard using historical IBTrACS tracks.        |
| IBTrACS_surge.ipynb         | Jupyter Notebook  | Example workflow for storm surge hazard using historical IBTrACS tracks.     |
| IBTrACS_wind.ipynb          | Jupyter Notebook  | Example workflow for wind hazard using historical IBTrACS tracks.            |
| STORM_surge.ipynb           | Jupyter Notebook  | Example workflow for storm surge hazard using synthetic STORM tracks.        |
| STORM_wind.ipynb            | Jupyter Notebook  | Example workflow for wind hazard using synthetic STORM tracks.               |


### Core Functions in TC_hazard_modelling.py
The main logic for the hazard computation is housed in the TC_hazard_modelling.py file.

**Centroids & Setup**
| Function                     | Description                                                                 |
|------------------------------|-----------------------------------------------------------------------------|
| create_centroids_from_shapefiles | Generates a list of CLIMADA Centroids (grid points) within the shapefile boundaries, applying an optional buffer and defined resolution. |
| plot_all_grids_with_zoom     | Visualizes the generated centroid grids and underlying shapefiles, with an optional zoom-in inset. |
| select_tracks                | Filters a TCTracks object to include only storms that pass through a specified geographic bounding box. |
| extract_selected_tracks      | Divides a multi-track TCTracks object into a dictionary of single-track TCTracks objects, keyed by Storm ID (SID). |


**Wind Hazard**
| Function          | Description                                                                 |
|-------------------|-----------------------------------------------------------------------------|
| compute_windfields | Computes the maximum wind speed at each centroid for all selected storm tracks using `climada.hazard.TropCyclone`. Saves the result to a CSV file. |


**Storm Surge Hazard**
| Function                  | Description                                                                 |
|---------------------------|-----------------------------------------------------------------------------|
| plot_windsurgerelation_Xu | Helper function to visualize the simplified wind-surge relationship used in the Bathtub model. |
| compute_surge             | Calculates the maximum storm surge height at each centroid using the `TCSurgeBathtub` model (a CLIMADA-Petals method), which requires an elevation file (`elev_file`). Saves the results to a CSV file. |


**Precipitation Hazard**
| Function                  | Description                                                                 |
|---------------------------|-----------------------------------------------------------------------------|
| downloading_ERA5_CDSAPI   | Downloads necessary ERA5 reanalysis data (Temperature, U/V Wind components at 600hPa/850hPa) for each storm from the Copernicus Climate Data Store (CDS). |
| extract_ERA5_parameters   | Processes the downloaded ERA5 data to extract the storm-specific T600, U850, and V850 parameters required for the TCRain model, and saves them to individual CSV files. |
| compute_precipitation     | Computes the total precipitation accumulation at each centroid using the `TCRain` model (CLIMADA-Petals method), which requires the pre-extracted ERA5 parameters attached to the track objects. |


### Getting Started
Dependencies: Ensure you have the required packages installed. The main dependencies include climada, climada_petals, numpy, pandas, geopandas, matplotlib, cartopy, and cdsapi.
Example installation command (may need adjustment based on your environment):
pip install numpy pandas geopandas matplotlib cartopy cdsapi climada climada_petals xarray


- **CDS API Key**: To download ERA5 data, you must configure your CDS API key as described in the Copernicus documentation.
- **Run Notebooks**: Start by exploring the five Jupyter notebooks (IBTrACS_*.ipynb and STORM_*.ipynb) to understand the step-by-step workflow for each hazard type. The notebooks demonstrate how to call the functions defined in TC_hazard_modelling.py. All notebooks follow the same structure, which is highlighted at the top of the notebook under _Overview_.
- **Data Placement**: Ensure your input files (shapefiles, tracks, DEMs) are correctly placed in the corresponding input folders (0-level-Leeward++, SRTM15Plus, STORM_PresentDayClimate) before running the notebooks.