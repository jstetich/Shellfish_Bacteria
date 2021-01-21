# Shellfish_Bacteria

Analysis of several years of data on bacteria levels in shellfish waters, from
Maine's Shellfish Sanitation Program

# Introduction
This archive includes R scripts and related products analyzing levels of
bacteria in waters where shellfish are harvested in Casco Bay, Maine.  The data
was provided  by the Department of Marine Resources, which is responsible for
tracking potential hazards to human health via consumption of shellfish.

MArine bivalve mollusks like clams, mussels, and oysters are filter feeders.
They ingest planktpn and bacterioplankton from their environment, and under
certain circumstances can concnntrate toxins or pathogens from the marine
environment.  Bivalves with elevated levels of toxins or pathogens can pose a
serious risk to human health.

The Bureau of Public Health oversees the application of the National Shellfish
Sanitation Program (NSSP) in Maine. The Department of Marine Resources (DMR)
collects related data that helps ensure safety of consumption of shellfish from
Maine waters.  Data collected includes data related to potential human health
hazards from consuming shellfish, such as the presence of phytotoxins in
shellfish tissue, presence of potentially toxic phytoplankton in marine waters,
and presence of elevated levels of indicator bacteria in local wasters.

Elevated levels of indicator bacteria do not necessarily indicate an iminant
threat to human health, they are an indicator of likely contamination of coastal
waters with fecal material, especially from warm blooded organisms. While risks
of exposure to (human) pathogens are highest if pollution comes from human
sources, any source of such contamination can be problematic.

# Statement of Purpose
CBEP is committed to the ideal of open science.  Our State of the Bay data
archives ensure the science underlying the 2020 State of the Bay report is
documented and reproducible by others. The purpose of these archives is to
release raw data and data analysis code whenever possible to allow others to
review, critique, learn from, and build upon CBEP science.

# Archive Structure
CBEP 2020 State of the Bay data analysis repositories are divided into from two
to four sub-folders.  All archives contain at least an "Original_Data" and a
"Graphics" folder.  The other two folders are only included if strictly
necessary.

- Original Data.  Original data, with a "DATA_SOURCES.md" or "READ ME.txt" 
file that documents data sources.
**DATA IN THIS FOLDER IS AS ORIGINALLY PROVIDED OR ACCESSED.** 

- Derived Data.  Data derived from the original raw data.  Includes
documentation of data reorganization steps, either in the form of files (R
notebooks, Excel files, etc.) that embody data transformations, or via README.md
or DATA_NOTES.md files.

- Analysis.  Contains one or more R Notebooks proceeding through the data
analysis steps. This often includes both preliminary data analysis --
principally graphical, and detailed analysis where necessary.

- Graphics.  Contains R Notebooks stepping through development of graphics, and
also copies of resulting graphics, usually in \*.png and \*.pdf formats.  These
graphics may differ from graphics as they appear in final State of the Bay
graphical layouts.

# Summary of Data Sources
We accessed publically available data on bacteria levels in Maine coastal waters
from the web site of Maine DMR.  These public data sets provide a synopsis of
recent observations of bacteria levels, usually from the most recent 30
observations collected from any particular location.  The summary data includes
median and 90th percentiles of those most recent observations.  These summary
values relate directly to Maine and NSSP safety standards.

Bacteria data is notoriously variable, with a handful of very high observations
scattered in within a much larger number of low values.  Thus certain
statistical analyses of bacteria levels requires access to raw data, not these
sumamries.  In late 2019, We requested data directly from maine DMR, and
recieved raw data  for the years 2015 through 2019 via e-mail.  All analyses are
based on these original observations.