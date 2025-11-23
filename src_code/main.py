def add2numbers(x, y):
    return x+y
print("KAshanHyderJaffri")
print("TAFENSW")
print( "KashanHyderJaffri")
print("TAFENSW")
print("Task1completed")
def myfunction (x,y,z):
    s=x+y+z
    return s
d = myfunction (1,2,3)


import math

def haversine_distance(lat1, lon1, lat2, lon2):
    # Radius of Earth in kilometers
    R = 6371.0

    # Convert degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    # Differences
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Haversine formula
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Distance in km
    distance = R * c
    return distance

# Sydney and Melbourne coordinates
lat_syd, lon_syd = 33.8688, 151.2093
lat_melb, lon_melb = 37.8136, 144.9631

# Calculate distance
distance_km = haversine_distance(lat_syd, lon_syd, lat_melb, lon_melb)

print("Haversine Distance between Sydney and Melbourne:", round(distance_km, 2), "km")

#import geopandas as gpd

#FILE_PATH = R"D:\Programming\ICTPRG440\spatial_data_original\groundwatertransitionalplans-shp"
#gdf = gpd.read_file(FILE_PATH)
#print (gdf.columns)



# Programmer: BBQ
# Creation Date: 19/10/2025
# Company: Gelos Enterprises
# Client: Equinox Surveyors
# Version number 3.9999
 
"""
Purpose of program.
    The client, Equinox Surveyors, has asked for a program to read and check the coordinate system of some given data.
    This data is then to be transformed to GDA 2020 ZOne 56 and exported in Esri shapefile format.
    The supplied information is provided in kml format. the client has specified that the reading and writing of this dta can be hardwired into the program.
    This program has gone through several iterations with more functionality added with each version.
    Password protection has been enabled
"""
'''
import os
 
def Passwordaccess():               #Creates password function.  
    Password = "555666"             #Sets password.
    guess = ""
    counter = 1                     #Sets guess counter to 1.
 
    while guess != Password and counter < 4:    # sets number of guesses to 3.
        guess = input(f"Enter Password, 3 tries allowed, Attempt {counter} - ")   #Text sent to command line, if entered password is incorrect it will loop 3 times.
        counter += 1
    if guess == Password:
        print("Access Granted")     #Displayed text When a correct password is given.
        return True
    else:
        print("You have exceeded the maximum allowed attempts! BYE BYE")            #Displayed text When 3 failed password attempts are made.
        return False
 
if __name__ == "__main__":                  #Protects main loop from running on import
 
    if  not Passwordaccess():                   #Invokes password. Modified from original code for looping.
        exit()
    while True:
'''
 
import geopandas as gpd         #GeoPandas library is used for reading and exporting geospatial vector data. gpd is alias.
 
'''
        The next lines define the location of the data file to be read, verified as kml and then location of the output file.
        this section has been modified from the original code from having the data file input and output hard coded to allowing the user to enter the file names.
        This allows for more lexibility and reuse of this program with different files.
        '''
   
 
filename = input("Enter the name of the file to import - include extension ")
input_dir = r"C:\TAFE Programming\ICTPRG440-Assignemt\Input"
INPUT_DATA = os.path.join(input_dir, filename)
           
if not os.path.exists(INPUT_DATA):                              #Checks if file exisist.
            print("Error: The input file does not exist.")
            exit()
 
if not INPUT_DATA.lower().endswith(".kml"):                     #Check if the input file is a KML file using operating system module.
                print("Error: The input file is not a KML file.","\n")       #If file is not a *.kml file then program terminates.
                exit()
else:
                print("data comfirmed as kml","\n")
                input("Press Enter to continue...",)
 
output_name = input("Enter the base name for the output shapefile (no extension): ")
output_dir = r"C:\TAFE Programming\ICTPRG440-Assignemt\Output"
output_filename = output_name + "_GDA2020_MGA_zone_56.shp"
OUTPUT_PATH = os.path.join(output_dir, output_filename)
 
if os.path.exists(OUTPUT_PATH):                                             #Places a warning if output already exists.
            print("Warning: Output file already exists and will be overwritten!")
            input("Press Enter to continue or Ctrl+C to cancel...")
   
 
gdf = gpd.read_file(INPUT_DATA) # Loads data into a GeoDataFrame (gdf) for processing.
 
 
def ShowAttributes(gdf): #Displays attribute table row by row
       
            print("Attribute Table:","\n")
            for index, row in gdf.iterrows():
                print(row)
 
ShowAttributes(gdf)
input("Press Enter to continue...")
 
print("Column headers:", gdf.columns.tolist())  #This lists the column headers from the kml file. Shape files, which will be the output have a 10 charater limit.
input("Press Enter to continue...")             #Pauses the program so the headers can be read.
 
EPSG = None                     # Initialise EPSG to None before querying original data file.
 
"""
        Next checks if the spatial data (gdf) has a defined coordinate reference system (CRS)
        Extracts EPSG code and prints to command line.
        If no CRS is found or the EPSG extraction fails, it prints a fail message.
        It then displays the first few rows of the data.
        """
       
if gdf.crs:
            try:
                EPSG = gdf.crs.to_epsg()
            except Exception:
                EPSG = None
 
if EPSG is not None:
            print("Original Coordinate system: EPSG:" + str(EPSG),"\n")
else:
            print("nil coordinate system present")
 
print(gdf.head(),"\n")  #Show first five rows of the gdf. head() is default 5 lines.
                                #The coordinates can be seen displayed here in Lat/long.
 
 
import warnings                                                 #Column names longer than 10 characters will be truncated, this supresses the error messages from this.
warnings.filterwarnings("ignore", category=UserWarning)         #Suppressing warnings is not ideal but the client insisted on using a shape file export.
warnings.filterwarnings("ignore", category=RuntimeWarning)      #This can be avoided by altering the original KML file.
 
def CordsystemTransformation(FilePath,EPSG=7856):               #Transforms the coordinate system
 
            gdf = gpd.read_file(FilePath)
            gdf = gdf.to_crs(epsg=EPSG)                                 #Reprojects the data to the specified EPSG coordinate system
            return gdf
 
Transformed = CordsystemTransformation(INPUT_DATA,EPSG=7856)    #Invokes the Coordinate transformation to GDA2020 / MGA zone 56
print(Transformed.head(),"\n")                                  #Shows first five rows of the transformed data. Coordinates can be seen and verified after transformation.
print("Coordinate System is now GDA2020 / MGA zone 56","\n")     #Prints Confirmation that the program has run and transformed the data.
Transformed.to_file(OUTPUT_PATH, driver="ESRI Shapefile")       #Creates ESRI Shapefile and saves in Output directory.
 
print("Output saved to:","\n", OUTPUT_PATH)                              #conformation of output location
 
 
again = input("Would you like to convert another file? (yes/no): ").strip().lower()     #Repeats if user wants to process another file.
if again != "yes":
            print("Have a nice day!","\n")
            #else:
            #break              




# main.py
"""
main.py
Author: Syed Kashan Hyder Jaffri (Jaffri)
Date: 2025-11-23
Purpose:
    - Read vector spatial data into a GeoDataFrame
    - Display the attribute table row-by-row in the console
    - Reproject vector spatial data from WGS84 to GDA2020 / MGA Zone 56 (EPSG:7856)
    - Save the reprojected layer as a shapefile in output/ directory
Notes:
    - Designed for GeoPandas 1.0.1
    - Test-mode: if input file not found, create a small sample GeoDataFrame to demonstrate functionality.
"""

import os
import sys
from pathlib import Path
import logging
import geopandas as gpd
from shapely.geometry import Point

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Constants
DEFAULT_OUTPUT_DIR = Path("output")
DEFAULT_OUTPUT_SHAPE = DEFAULT_OUTPUT_DIR / "projected_layer.shp"
# EPSG codes: input WGS84 (EPSG:4326) -> GDA2020 / MGA Zone 56 (EPSG:7856)
TARGET_CRS_EPSG = 7856
SOURCE_CRS_EPSG = 4326


def read_vector_to_gdf(filepath):
    """
    Read a vector spatial file and return a GeoDataFrame.

    Parameters:
        filepath (str or pathlib.Path) : Path to the vector file (shapefile, geojson, etc.)
            - Data type: str or pathlib.Path

    Returns:
        gdf (geopandas.GeoDataFrame) : GeoDataFrame containing the input layer
            - Data type: geopandas.GeoDataFrame
        If an error occurs, returns None.

    Description:
        Uses geopandas.read_file to load vector spatial data. Wraps call in try/except
        to prevent the program from crashing if the file doesn't exist or is unreadable.
    """
    try:
        filepath = Path(filepath)
        logging.info(f"Attempting to read vector file: {filepath}")
        if not filepath.exists():
            logging.warning(f"File not found: {filepath}")
            return None
        gdf = gpd.read_file(str(filepath))
        logging.info(f"Loaded layer with {len(gdf)} features and CRS: {gdf.crs}")
        return gdf
    except Exception as e:
        logging.error(f"Error reading vector file: {e}")
        return None


def print_attributes_row_by_row(gdf):
    """
    Print the attribute table row-by-row to the console.

    Parameters:
        gdf (geopandas.GeoDataFrame) : The GeoDataFrame whose rows will be printed
            - Data type: geopandas.GeoDataFrame

    Returns:
        None

    Description:
        Iterates through the rows of the GeoDataFrame using a loop (for ... in).
        For each row, prints the index and a dictionary of attribute values (excluding geometry).
        The function uses try/except to handle unexpected errors.
    """
    try:
        if gdf is None:
            logging.warning("No GeoDataFrame provided to print.")
            return None

        if gdf.empty:
            logging.info("GeoDataFrame is empty; nothing to print.")
            return None

        # Use a for loop to iterate row-by-row
        for idx, row in gdf.iterrows():
            # Exclude the geometry column from attribute printout for readability
            attributes = {col: row[col] for col in gdf.columns if col != gdf.geometry.name}
            logging.info(f"Row {idx}: {attributes}")
        return None
    except Exception as e:
        logging.error(f"Error printing attributes: {e}")
        return None


def reproject_and_save(gdf, target_epsg=TARGET_CRS_EPSG, out_shapefile=DEFAULT_OUTPUT_SHAPE):
    """
    Reproject a GeoDataFrame from its source CRS (or WGS84 if undefined) to a target CRS,
    and save the reprojected layer as a shapefile.

    Parameters:
        gdf (geopandas.GeoDataFrame) : Input GeoDataFrame to reproject.
            - Data type: geopandas.GeoDataFrame
        target_epsg (int) : EPSG code of the desired projected coordinate system.
            - Data type: int
            - Default: TARGET_CRS_EPSG (7856 for GDA2020 / MGA zone 56)
        out_shapefile (str or pathlib.Path) : Path to the output shapefile to write.
            - Data type: str or pathlib.Path
            - Default: output/projected_layer.shp

    Returns:
        out_path (pathlib.Path) : Path object pointing to the saved shapefile if successful.
            - Data type: pathlib.Path
        If an error occurs, returns None.

    Description:
        Checks the input GeoDataFrame for a CRS. If it's None, assumes WGS84 (EPSG:4326).
        Then calls GeoDataFrame.to_crs to reproject and saves using GeoDataFrame.to_file.
        The function ensures the output folder exists, uses try/except to catch errors,
        and includes an 'if' check to avoid reprojecting an empty GeoDataFrame.
    """
    try:
        if gdf is None:
            logging.error("No GeoDataFrame provided for reprojection.")
            return None

        if gdf.empty:
            logging.warning("GeoDataFrame is empty. Nothing to reproject.")
            return None

        # Ensure output directory exists
        out_shapefile = Path(out_shapefile)
        out_dir = out_shapefile.parent
        out_dir.mkdir(parents=True, exist_ok=True)

        # Ensure the gdf has a CRS; if not, assume WGS84 (EPSG:4326)
        if gdf.crs is None:
            logging.warning("Input GeoDataFrame has no CRS. Assuming WGS84 (EPSG:4326).")
            gdf = gdf.set_crs(epsg=SOURCE_CRS_EPSG, allow_override=True)

        logging.info(f"Reprojecting from {gdf.crs} to EPSG:{target_epsg} ...")
        gdf_proj = gdf.to_crs(epsg=target_epsg)

        # Save to shapefile
        logging.info(f"Saving reprojected layer to {out_shapefile} ...")
        gdf_proj.to_file(str(out_shapefile))
        logging.info("Save completed.")
        return out_shapefile
    except Exception as e:
        logging.error(f"Error reprojecting/saving layer: {e}")
        return None


def create_sample_gdf():
    """
    Create a small sample GeoDataFrame (point features) in WGS84 for testing purposes.

    Parameters:
        None

    Returns:
        sample_gdf (geopandas.GeoDataFrame) : A GeoDataFrame with two sample points (Melbourne, Sydney)
            - Data type: geopandas.GeoDataFrame

    Description:
        Convenience helper to create a test layer when no input file is provided.
        The CRS will be set to EPSG:4326 (WGS84).
    """
    try:
        # Sydney and Melbourne points (lon, lat)
        data = {
            "name": ["Sydney", "Melbourne"],
            "latitude": [33.8688, 37.8136],
            "longitude": [151.2093, 144.9631],
        }
        geometries = [Point(xy) for xy in zip(data["longitude"], data["latitude"])]
        sample_gdf = gpd.GeoDataFrame(data, geometry=geometries, crs=f"EPSG:{SOURCE_CRS_EPSG}")
        return sample_gdf
    except Exception as e:
        logging.error(f"Error creating sample GeoDataFrame: {e}")
        return None


def main(input_path=None):
    """
    Main driver function to run read -> print -> reproject -> save workflow.

    Parameters:
        input_path (str or pathlib.Path, optional) : Path to the vector input file.
            - Data type: str or pathlib.Path
            - If None or file missing, a sample GeoDataFrame will be used.

    Returns:
        None

    Description:
        Executes the workflow. Uses try/except to handle top-level errors and returns exit codes
        to indicate success (0) or failure (1). Demonstrates use of sequence, selection, and iteration.
    """
    try:
        # Attempt to read the provided vector file
        gdf = None
        if input_path:
            gdf = read_vector_to_gdf(input_path)

        # If read failed or no input provided, create sample data for demonstration/testing
        if gdf is None:
            logging.info("Using sample GeoDataFrame for testing.")
            gdf = create_sample_gdf()

        # Print attributes row by row (loop requirement)
        print_attributes_row_by_row(gdf)

        # Reproject to target CRS and save
        out_file = reproject_and_save(gdf, target_epsg=TARGET_CRS_EPSG)
        if out_file:
            logging.info(f"Projected shapefile saved at: {out_file}")
        else:
            logging.error("Failed to save projected shapefile.")
        return 0
    except Exception as e:
        logging.critical(f"Unhandled exception in main: {e}")
        return 1


if __name__ == "__main__":
    # Optionally accept a file path argument from command line
    input_fp = None
    if len(sys.argv) > 1:
        input_fp = sys.argv[1]

    exit_code = main(input_fp)
    sys.exit(exit_code)