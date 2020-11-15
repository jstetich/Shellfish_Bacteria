# Data Notes For Shellfish Sanitation Data
## Data Attributes from "P90" files
We lack complete metadata for these files,a nd so have to interpret some fields.

After merging the three p90 data files and cleaning up inconsistencies, we found
the following data columns that were consistent across all three years. Our
interpretation of data columns is provided in in parentheses.

*  "Station" (Station Code -- related to GIS data by lat long, and for 2018
   data, by x and y in UTM)
*  "Class" (Approved, Conditionally Approved, Restricted, Conditionally
   Restricted or Prohibited).  The official harvesting status of the area where
   the sampling site is located.
*  "Count" (Number of Samples. ordinarily 30, sometimes less, because DMR
   calculates the "median" and "p90" values based on the most recent 30 samples.
*  "MFCount" (Some times lower that Count. Why? Which number actually represents
   the number of samples included in calculations?)
*  "GM" (Geometric Mean)
*  "SDV" (Some sort of a standard deviation? Is this SD of raw (untransformed)
   data, SD of log data, or something else?) 
*  "MAX" (Maximum Observed Value)
*  "P90" (90th percentile)
*  "Appd_Std" (Standard for a station to be in "approved" condition. 
   Almost always 31)
*  "Restr_Std" (Standard for a station to be in "restricted" condition.
   Almost always 163)
*  "Lat_DD" (Latitude. Assumed to be in WGS84)
*  "Long_DD" (Longitude. Assumed to be in WGS84)
*  "Grow_Area" (Code for management areas for classification and closure
   decisions).


## GIS Data
GIS data were assembled from the p90 data, which included lats and longs.
Statewide data was read into ArcGIS as a CSV table, and saved as a shapefile.
Locations were selected that were within 500 m of the Casco Bay waterhsed layer,
which go al lsampling locations in our region.

Reviewing the data, however, amkes it prety clear the p90 data in subsequent years are not independent, so other than using the 2018 p90 ma, we wil lneed toanalyze underlying data to get meaningful trend analyses.