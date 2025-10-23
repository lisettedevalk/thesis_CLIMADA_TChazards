# ---- INSTALL REQUIRED PACKAGES ----

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from shapely.geometry import box
from shapely.geometry import Point
import cdsapi

from climada.hazard import TropCyclone, TCTracks, Centroids
from climada_petals.hazard import TCSurgeBathtub
from climada_petals.hazard import TCRain



# ---- CREATE CENTROIDS ----

def create_centroids_from_shapefiles(folder_path, resolution_km, buffer_km):
    """
    Create a list of CLIMADA Centroids objects from shapefiles in a folder,
    with an added buffer around each shape.

    Parameters
    ----------
    folder_path : str
        Path to folder containing .shp files.
    resolution_km : float
        Spacing between grid points in km.
    buffer_km : float
        Distance in km to buffer outward from the polygon boundary.

    Returns
    -------
    centroids_list : list
        List of CLIMADA Centroids objects, one per shapefile.
    """
    centroids_list = []
    
    # Store shapefile names for later use
    shapefiles_list = []

    for filename in os.listdir(folder_path):
        if filename.endswith(".shp"):
            filepath = os.path.join(folder_path, filename)
            gdf = gpd.read_file(filepath).to_crs("EPSG:4326")

            # Reproject to meters, buffer in km, then back to degrees
            gdf_m = gdf.to_crs(epsg=3857)  # meters
            gdf_buffered = gdf_m.buffer(buffer_km * 1000).to_crs("EPSG:4326")

            # Merge shapes into a single polygon
            gdf_union = gdf_buffered.unary_union
            
            # Generate points within the buffered polygon
            minx, miny, maxx, maxy = gdf_union.bounds
            
            # Create a grid of points
            step_deg = resolution_km / 111.0 # ~1 degree ≈ 111 km near equator
            x_coords = np.arange(minx, maxx, step_deg)
            y_coords = np.arange(miny, maxy, step_deg)
            
            grid_points = []
            for x in x_coords:
                for y in y_coords:
                    grid_points.append(Point(x, y))

            grid_gdf = gpd.GeoDataFrame(geometry=grid_points, crs="EPSG:4326")

            # Select points that fall within the buffered polygon
            points_in_polygon = grid_gdf[grid_gdf.within(gdf_union)]

            if points_in_polygon.empty:
                print(f"Warning: No grid points found for {filename} after buffering.")
                continue

            # Convert to CLIMADA Centroids
            centroids_list.append(Centroids.from_geodataframe(points_in_polygon))
            shapefiles_list.append(filename) # Store the filename

    print(f"Created centroids for {len(centroids_list)} shapefiles.")

    return centroids_list, shapefiles_list



def create_coastal_centroids_from_shapefiles(folder_path, resolution_km, band_width_km, inward_buffer_km=0.5):
    """
    Create a list of CLIMADA Centroids objects for a coastal band
    surrounding shapefiles, excluding the interior of the country shape.

    Parameters
    ----------
    folder_path : str
        Path to folder containing .shp files.
    resolution_km : float
        Spacing between grid points in km.
    band_width_km : float
        Width of the coastal band in km (e.g., 5.0 km, which is mostly offshore).
    inward_buffer_km : float, optional
        Distance in km to buffer inward from the coastline (e.g., 0.5 km).

    Returns
    -------
    centroids_list : list
        List of CLIMADA Centroids objects, one per shapefile.
    shapefiles_list : list
        List of shapefile names.
    """
    centroids_list = []
    shapefiles_list = []
    
    outer_buffer_km = band_width_km + inward_buffer_km 

    for filename in os.listdir(folder_path):
        if filename.endswith(".shp"):
            filepath = os.path.join(folder_path, filename)
            gdf = gpd.read_file(filepath).to_crs("EPSG:4326")
            
            # Reproject to a projected CRS (meters) for accurate buffering
            gdf_m = gdf.to_crs(epsg=3857)  # meters

            
            # This polygon represents the country + the entire coastal band.
            outer_buffer_m = gdf_m.buffer(outer_buffer_km * 1000).unary_union
            
            # This polygon represents the country minus the coastal strip on land.
            inner_buffer_m = gdf_m.buffer(-inward_buffer_km * 1000).unary_union
            
            # This removes the interior of the country, leaving the ring shape.
            gdf_band = outer_buffer_m.difference(inner_buffer_m)
            
            # Reproject the resulting band back to degrees for grid generation
            gdf_band_gpd = gpd.GeoSeries(gdf_band, crs=3857).to_crs("EPSG:4326")
            gdf_union = gdf_band_gpd.unary_union # Ensure it's treated as a single geometry

            
            minx, miny, maxx, maxy = gdf_union.bounds
            
            step_deg = resolution_km / 111.0 
            x_coords = np.arange(minx, maxx, step_deg)
            y_coords = np.arange(miny, maxy, step_deg)
            
            grid_points = []
            for x in x_coords:
                for y in y_coords:
                    grid_points.append(Point(x, y))

            grid_gdf = gpd.GeoDataFrame(geometry=grid_points, crs="EPSG:4326")

            # Select points that fall within the coastal band (gdf_union)
            points_in_polygon = grid_gdf[grid_gdf.within(gdf_union)]

            if points_in_polygon.empty:
                print(f"Warning: No grid points found for {filename} after coastal band creation.")
                continue

            centroids_list.append(Centroids.from_geodataframe(points_in_polygon))
            shapefiles_list.append(filename) 

    print(f"Created centroids for {len(centroids_list)} shapefiles.")

    return centroids_list, shapefiles_list






def plot_all_grids_with_zoom(folder_path, shapefiles_list, centroids_list, my_colors, resolution_km, zoom_extent):
    """
    Plot all shapefiles with their precomputed centroids on a single map,
    including a zoomed-in inset.

    Parameters
    ----------
    folder : str
        Path to the folder containing shapefiles.
    shapefiles_list : list
        List of shapefile filenames (must match centroids_list order).
    centroids_list : list
        List of CLIMADA Centroids objects (output of create_centroids_from_shapefiles).
    my_colors : dict
        Dictionary with color definitions, e.g., {"510darkblue": "#...", "510purple": "#..."}.
    resolution_km : float
        Grid resolution in kilometers (for displaying in plot title).
    zoom_extent : tuple, optional
        (xmin, xmax, ymin, ymax) for the zoomed-in inset. Default is focused near Anguilla.

    Returns
    -------
    None
        Displays a matplotlib figure.
    """
    if len(shapefiles_list) != len(centroids_list):
        print("Warning: shapefiles_list and centroids_list lengths differ. Skipping missing ones.")

    all_shapes = []
    all_points = []

    for filename, centroids in zip(shapefiles_list, centroids_list):
        if centroids.gdf.empty:
            continue  # skip empty centroid sets

        filepath = os.path.join(folder_path, filename)
        gdf = gpd.read_file(filepath).to_crs("EPSG:4326")
        iso3 = filename.split("_")[1]  # Extract ISO3 code (e.g., gadm41_AIA_0.shp → "AIA")
        gdf["ISO3"] = iso3

        all_shapes.append(gdf)
        all_points.append(centroids.gdf)

    # Combine all into single GeoDataFrames
    shapes_gdf = gpd.GeoDataFrame(pd.concat(all_shapes, ignore_index=True), crs="EPSG:4326")
    points_gdf = gpd.GeoDataFrame(pd.concat(all_points, ignore_index=True), crs="EPSG:4326")

    # --- Main plot ---
    fig, ax = plt.subplots(figsize=(18, 12))
    shapes_gdf.plot(ax=ax, edgecolor=my_colors["510darkblue"], facecolor="white")
    points_gdf.plot(ax=ax, markersize=0.5, color=my_colors["510purple"], alpha=0.7)

    # Add country labels
    offset_x = 0.08
    offset_y = 0.08
    for _, row in shapes_gdf.iterrows():
        centroid = row.geometry.centroid
        ax.text(
            centroid.x + offset_x,
            centroid.y + offset_y,
            row["ISO3"],
            ha="center", va="center",
            fontsize=14, fontweight="bold",
            color="black",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.5, pad=1)
        )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Longitude (°)", fontsize = 16)
    ax.set_ylabel("Latitude (°)", fontsize = 16)
    ax.set_title(f"A {resolution_km} km grid over all islands", fontsize=18)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 14)

    # --- Draw zoom box ---
    zoomxmin, zoomxmax, zoomymin, zoomymax = zoom_extent
    ax.plot(
        [zoomxmin, zoomxmax, zoomxmax, zoomxmin, zoomxmin],
        [zoomymin, zoomymin, zoomymax, zoomymax, zoomymin],
        color="black", linewidth=1, linestyle="--"
    )

    # --- Create zoomed inset ---
    axins = inset_axes(
        ax, width="60%", height="60%", loc="lower left",
        bbox_to_anchor=(-0.06, 0.05, 1, 1),
        bbox_transform=ax.transAxes, borderpad=1
    )
    shapes_gdf.plot(ax=axins, edgecolor=my_colors["510darkblue"], facecolor="white")
    points_gdf.plot(ax=axins, markersize=5, color=my_colors["510purple"], alpha=0.8)

    # Add country labels in zoom view (only those inside box)
    for _, row in shapes_gdf.iterrows():
        centroid = row.geometry.centroid
        if zoomxmin <= centroid.x <= zoomxmax and zoomymin <= centroid.y <= zoomymax:
            axins.text(
                centroid.x + offset_x,
                centroid.y + offset_y / 2,
                row["ISO3"],
                ha="center", va="center",
                fontsize=14, fontweight="bold",
                color="black",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.5, pad=1)
            )

    axins.set_xlim(zoomxmin, zoomxmax)
    axins.set_ylim(zoomymin, zoomymax)
    axins.set_aspect("equal", adjustable="box")
    ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
    


# ---- PRE-PROCESS TRACKS ----

def select_tracks(tracks, lon_min, lon_max, lat_min, lat_max):
    """
    Select tropical cyclone tracks that pass through a given bounding box
    and interpolate them to a specified time step.

    Parameters
    ----------
    tracks : TCTracks
        A TCTracks object containing the original storm data.
    lon_min, lon_max, lat_min, lat_max : float
        Bounding box coordinates (in degrees).

    Returns
    -------
    TCTracks
        A new TCTracks object with storms inside the bounding box,
        interpolated to the specified time step.
    """
    selected_tracks = []
    
    for track in tracks.data:
        lons = np.array(track.lon)
        lats = np.array(track.lat)

        # Check if storm enters the box
        if np.any((lon_min <= lons) & (lons <= lon_max) &
                  (lat_min <= lats) & (lats <= lat_max)):
            selected_tracks.append(track)

    # Create new TCTracks object
    filtered_tracks = TCTracks()
    filtered_tracks.data = selected_tracks

    return filtered_tracks


def extract_selected_tracks(tracks):
    """
    Create a dictionary of TCTracks objects, one per unique SID.

    Parameters
    ----------
    tracks : TropCyclone
        A TropCyclone object with .data containing xarray track datasets.
        
    Returns
    -------
    dict
        Dictionary with SID as key and TCTracks (with single track) as value.
    """
    # Collect all SIDs
    sid_list = [track.attrs.get('sid') for track in tracks.data]
    print("Number of tracks:", len(set(sid_list)))

    # Build dictionary of TCTracks per SID
    tr_dict = {}
    for track in tracks.data:
        sid = track.attrs['sid']  # or track.sid in older versions
        new_tc = TCTracks()
        new_tc.data = [track]  # keep only this track
        tr_dict[sid] = new_tc

    return tr_dict




# ---- HAZARD: WIND ----

def compute_windfields(tc_tracks_dict, centroids_list, shapefiles, output_csv_path):
    """
    Compute wind fields for multiple TropCyclone objects and save results.

    Parameters:
    - tc_tracks_dict: dict with keys=sid, values=TCTracks objects
    - centroids_list: list of Centroids objects
    - shapefiles: list of filenames like ['grid_ATG.shp', ...]
    - output_csv_path: str path to output CSV
    """

    records = []

    for sid, tr in tc_tracks_dict.items():
        for i, filename in enumerate(shapefiles):
            iso = filename.split('_')[1].split('.')[0]  # Remove .shp extension
            try:
                centroids = centroids_list[i]
                tc_result = TropCyclone.from_tracks(tr, centroids, intensity_thres = 17.5)

                zvals = tc_result.intensity.toarray()[0]
                for lon, lat, wind in zip(centroids.lon, centroids.lat, zvals):
                    records.append({
                        'SID': sid,
                        'ISO': iso,
                        'lon': lon,
                        'lat': lat,
                        'max_wind_speed_mps': wind
                    })

            except Exception as e:
                print(f"Error for SID {sid} and ISO {iso}: {e}")
                continue

    df = pd.DataFrame(records)
    df.set_index(['SID', 'ISO'], inplace=True)
    df.to_csv(output_csv_path)
    print(f"All results saved to {output_csv_path}")
    return df



# ---- HAZARD: SURGE ----

def plot_windsurgerelation_Xu(mph2ms = 0.44704, f2m = 0.3048):
    # the points read from the SLOSH graph
    v0 = 60*mph2ms;
    v1 = 140*mph2ms;
    s0 = 6*f2m;
    s1 = 18*f2m;
    
    # the parameters for the linear function: a*(v-v0)+s0
    a = (s1-s0)/(v1-v0)

    # graphical representation
    v = np.arange(20, 100)
    vmph = np.arange(60, 141, 20)
    
    plt.plot(v, a*(v-v0)+s0, label='approximation')
    plt.scatter(vmph*mph2ms, a*(vmph*mph2ms-v0)+s0, label='SLOSH points')
    plt.xlabel('Wind speed (m/s)')
    plt.ylabel('Surge height (m)')
    plt.title('Wind-surge relationship in Xu (2010)')
    plt.grid()
    plt.legend()




def compute_surge(tc_tracks_dict, centroids_list, shapefiles, elev_file, output_csv_path, intensity_thres=17.5):
    """
    Calculate storm surge for multiple storms (all countries together) and save in one file.

    Parameters
    ----------
    tc_tracks_dict : dict
        Dictionary with keys = SID, values = TCTracks objects.
    centroids_list : list
        List of Centroids objects (one per country).
    shapefiles : list
        List of shapefile paths corresponding to centroids.
    elev_file : str
        Path to the elevation file.
    intensity_thres : float, optional
        Wind speed threshold (m/s) for selecting tracks (default 17.5).
    output_csv_path : str, optional
        Path where the combined CSV will be saved.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with surge results for all countries.
    """

    records = []

    for sid, tr in tc_tracks_dict.items():
        for i, filename in enumerate(shapefiles):
            iso = filename.split('_')[1].split('.')[0]  # Extract ISO
            try:
                centroids = centroids_list[i]
                tc = TropCyclone.from_tracks(tr, centroids, intensity_thres=intensity_thres)
                if tc.intensity.shape[0] == 0:
                    continue  # no storms above threshold

                tc_surge = TCSurgeBathtub.from_tc_winds(tc, elev_file)
                zvals = np.array(tc_surge.intensity.todense())[0].flatten()

                for lon, lat, surge in zip(centroids.lon, centroids.lat, zvals):
                    records.append({
                        'SID': sid,
                        'ISO': iso,
                        'lon': lon,
                        'lat': lat,
                        'surge_m': surge
                    })
                
            except Exception as e:
                print(f"Error for SID {sid} and ISO {iso}: {e}")
                continue

    #######################

    with open(f"{output_csv_path}.txt", "w") as f:
        for item in records:
            f.write(f"{item['SID'].split('_')[-1]},{item['ISO']},{item['lon']:.5f},{item['lat']:.5f},{item['surge_m']:.2f}\n")

    #######################
    
    ## build final DataFrame
    #df = pd.DataFrame(records)
    #
    ## save one combined CSV
    #df.to_csv(output_csv_path, index=False)
    #print(f"Combined surge results saved to {output_csv_path}")
    #
    #return df





# ---- HAZARD: PRECIPITATION ----

def downloading_ERA5_CDSAPI(sid_list, tc_track_dict, lon_min, lon_max, lat_min, lat_max, output_nc_dir):
    """
    Download ERA5 pressure level data for given storm IDs.

    Parameters
    ----------
    sid_list : list
        List of storm IDs.
    tr_ : dict
        Dictionary with track data (per storm ID).
    lon_min, lon_max, lat_min, lat_max : float
        Bounding box coordinates.
    output_nc_path : str
        Directory to save the NetCDF files.
    """

    # Initialize CDS API client once
    c = cdsapi.Client()

    for sid in sid_list:
        # Create date and time strings from storm track data
        date_strings = np.unique(tc_track_dict[sid].data[0].time.dt.strftime("%Y-%m-%d"))
        time_steps = np.unique(tc_track_dict[sid].data[0].time.dt.strftime("%H:%M"))

        # Define bounding box: [North, West, South, East]
        bounding_box = [lat_max, lon_min, lat_min, lon_max]

        # Create unique output filename
        output_csv_path = f'{output_nc_dir}/era5_data_{sid}.nc' 

        # Request parameters
        request_params = {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': [
                'temperature',
                'u_component_of_wind',
                'v_component_of_wind'
            ],
            'pressure_level': ['600', '850'],
            'date': list(date_strings),
            'time': list(time_steps),
            'area': bounding_box,
        }

        # Try to download, skip if fails
        try:
            c.retrieve('reanalysis-era5-pressure-levels', request_params, output_csv_path)
            print(f"Download complete for {sid}. Saved to {output_csv_path}")
        except Exception as e: 
            print(f"Skipping {sid} due to request failure: {e}")



def extract_ERA5_parameters(sid_list, tc_track_dict, input_nc_dir, output_csv_dir):
    """
    Extract ERA5 parameters (T600, U850, V850) for each storm in sid_list
    and save them as CSV files.

    Parameters
    ----------
    sid_list : list
        List of storm IDs.
    tr_ : dict
        Dictionary with storm track data.
    input_nc_path : str
        Directory containing ERA5 NetCDF files.
    output_csv_path : str
        Directory to save the parameter CSV files.
    """

    for sid in sid_list:
        # --- Create in- and output paths ---
        input_nc_path = f'{input_nc_dir}/era5_data_{sid}.nc'
        output_csv_path = f'{output_csv_dir}/ERA5_parameters_{sid}.csv'
        
        # --- Load NetCDF file ---
        ds = xr.open_dataset(input_nc_path)

        # --- Extract storm info ---
        storm_centre_lon = tc_track_dict[sid].data[0].lon
        storm_centre_lat = tc_track_dict[sid].data[0].lat
        time = tc_track_dict[sid].data[0].time

        # --- T600 ---
        t600 = ds["t"].sel(pressure_level=600)
        T600 = t600.sel(
            latitude=storm_centre_lat,
            longitude=storm_centre_lon,
            valid_time=time,
            method="nearest"
        )

        # --- U850 & V850 ---
        u850 = ds["u"].sel(pressure_level=850, valid_time=time, method="nearest")
        v850 = ds["v"].sel(pressure_level=850, valid_time=time, method="nearest")

        # Convert to radians
        lat_rad = np.deg2rad(u850.latitude)
        lon_rad = np.deg2rad(u850.longitude)
        storm_lat_rad = np.deg2rad(storm_centre_lat)
        storm_lon_rad = np.deg2rad(storm_centre_lon)

        # --- Haversine distance (km) ---
        R = 6371.0
        dlat = lat_rad - storm_lat_rad
        dlon = storm_lon_rad - lon_rad
        a = np.sin(dlat / 2) ** 2 + np.cos(storm_lat_rad) * np.cos(lat_rad) * np.sin(dlon / 2) ** 2
        c = 2 * np.arcsin(np.sqrt(a))
        dist_km = R * c

        # --- Apply annulus mask (200–500 km) ---
        mask = (dist_km >= 200) & (dist_km <= 500)
        dist_km_reordered = dist_km.transpose("time", "latitude", "longitude")
        mask_reordered = mask.transpose("time", "latitude", "longitude")

        u850_annulus = u850.where(mask_reordered)
        v850_annulus = v850.where(mask_reordered)

        # --- Compute annulus mean ---
        U850 = u850_annulus.mean(dim=("latitude", "longitude"))
        V850 = v850_annulus.mean(dim=("latitude", "longitude"))

        # --- Save to CSV ---
        era5_df = pd.DataFrame({
            "T600": T600.values,
            "U850": U850.values,
            "V850": V850.values
        }, index=pd.to_datetime(time.values))

        era5_df.to_csv(output_csv_path)
    print(f"All results saved, in {len(sid_list)} files.")


def compute_precipitation(sid_list, tc_track_dict, shapefiles_list, centroids_list, input_csv_path, model_kwargs, output_csv_path):
    """
    Compute TC rainfall using TCRain with ERA5 parameters attached to storm tracks.

    Parameters
    ----------
    sid_list : list
        List of storm IDs.
    tc_track_dict : dict
        Dictionary with storm track data (modified to accept ERA5 parameters).
    shapefiles_list : list
        List of shapefile names (used to extract ISO codes).
    centroids_list : list
        List of centroid objects corresponding to shapefiles.
    input_csv_path : str
        Directory where storm-specific ERA5 parameter CSVs are stored. (e.g., 'ERA5_parameters/')
    model_kwargs : dict, optional
        Dictionary of keyword arguments passed to the TCRain model.
    output_csv_path : str
        Path to save the final precipitation CSV.
    """

    records = []

    for sid in sid_list:
        # --- CORRECTED LINE: Load the storm-specific CSV ---
        # The storm-specific file is constructed by joining the base path and the SID.
        era5_file_path = os.path.join(input_csv_path, f"ERA5_parameters_{sid}.csv")
        
        try:
            era5_parameters = pd.read_csv(era5_file_path)
        except FileNotFoundError:
            print(f"Skipping SID {sid}: ERA5 parameter file not found at {era5_file_path}")
            continue

        tc_track_dict[sid].data[0]["t600"] = era5_parameters["T600"].values
        tc_track_dict[sid].data[0]["u850"] = era5_parameters["U850"].values
        tc_track_dict[sid].data[0]["v850"] = era5_parameters["V850"].values
        
        for i, filename in enumerate(shapefiles_list):
            iso = filename.split('_')[1].split('.')[0]
            try:
                centroids = centroids_list[i]
                TCR_ = TCRain.from_tracks(tc_track_dict[sid], 
                                          centroids, 
                                          model="TCR", 
                                          model_kwargs=model_kwargs,
                                          max_dist_eye_km = 800)
                zvals = np.array(TCR_.intensity.todense())[0].flatten()
        
                for lon, lat, precipitation in zip(centroids.lon, centroids.lat, zvals):
                    records.append({
                        'SID': sid,
                        'ISO': iso,
                        'lon': lon,
                        'lat': lat,
                        'prec_total_mm': precipitation
                    })
            except Exception as e:
                # The original error (size mismatch) occurred here, but with the fix 
                # it's now less likely and will catch true TCRain errors.
                print(f"Error for SID {sid} and ISO {iso}: {e}")
    
    df = pd.DataFrame(records)
    df.set_index(['SID', 'ISO'], inplace=True)
    df.to_csv(output_csv_path)
    print("All results saved.")