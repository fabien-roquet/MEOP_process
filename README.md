# MEOP_process
Scripts used to process MEOP data (meop.net)

Since 2004, several hundred thousands profiles of temperature and salinity have been 
collected by instrumented animals. The use of elephant seals has been particularly 
effective to sample the Southern Ocean and the North Pacific. These hydrographic data 
have been assembled in a quality-controlled database, the MEOP-CTD database, that can 
be accessed through this website.
For more information, visit the website meop.net
For any questions, contact info@meop.net or fabien.roquet@gmail.com


## THE MEOP-CTD DATABASE
README: the MEOP-CTD database (owner: Fabien Roquet and the MEOP consortium)
Release Date: 11/11/2017
Version name: MEOP-CTD_2017-11-11
https://opendatacommons.org/licenses/odbl/


## DATA FORMATS

Data are provided in three different formats. 
* _DATA_ncARGO:_ For a thourough scientific use of the data, or for oceanographic data centers, it is advised to use the 
marine mammal netCDF format (files in DATA_ncARGO) as it serves as the reference. This format can be 
easily read in Ocean Data View, using the Import/ARGO profiles/Float profiles menu, or using your 
favorite data processing software (e.g. Python, Matlab, IDL). 
* _DATA_ncARGO_interp:_ For ease of use, the DATA_ncARGO_interp provides the same data as in DATA_ncARGO, except it has
been interpolated on a regular vertical grid (1dbar spacing).
* _DATA_csv_interp:_ A csv format (ASCII) is also provided (files in DATA_csv_interp) which can be opened with Excel
or any text editor. Here, only data flagged as good are included, and are given on a regular 
vertical grid (1dbar spacing).



## DATA INFORMATION

The data that is publicly available is shown in the figure map_global_public.png. More data 
is available upon request, as part of the MEOP-CTD database. See map_global_private.png for
the distribution of private data.
Important metadata and statistics are listed in the info_*.csv files:
* _info_total.csv_ gives global statistics about the MEOP-CTD database
* _info_groups.csv_ gives statistics by national groups (see [MEOP groups](meop.net/groups/) for information)
* _info_deployments.csv_ gives statistics by deployment
* _info_tags.csv_ gives information and statistics by individual tag.
For each deployments, distribution maps are available in the MAPS directory, and a pdf
document with basic plots of CTD data (TS plots, time sections) is provided in the
directory PDF.



## HOW TO CITE

If you use this dataset for a publication, please add the following sentence 
in the Acknowledgement part:
"The marine mammal data were collected and made freely available by the International MEOP 
Consortium and the national programs that contribute to it (http://www.meop.net)."

Also consider citing the following papers when you use the MEOP-CTD dataset
for oceanographic applications:
- Treasure, A. M., Roquet, F., Ansorge, I. J., Bester, M. N., Horst Bornemann, L. B., Charrassin, J.-B., Chevallier, D., Costa, D. P., Fedak, M. A., Guinet, C., Hammill, M. O., Harcourt, R. G., Hindell, M. A., Kovacs, K. M., Lea, M.-A., Lovell, P., Lowther, A. D., Lydersen, C., McIntyre, T., McMahon, C. R., Muelbert, M. M. C., Nicholls, K., Picard, B., Reverdin, G., Trites, A. W., Williams, G. D., and de Bruyn, P. J. N., 2017. Marine Mammals Exploring the Oceans Pole to Pole: A Review of the MEOP Consortium. Oceanography, 30(2):132–138, doi: 10.5670/oceanog.2017.234
- Roquet F., Williams G., Hindell M. A., Harcourt R., McMahon C. R., Guinet C., Charrassin 
J.-B., Reverdin G., Boehme L., Lovell P. and Fedak M. A., 2014. A Southern Indian Ocean 
database of hydrographic profiles obtained with instrumented elephant seals. Nature 
Scientific Data, 1:140028, doi: 10.1038/sdata.2014.28

### Important technical papers : 
* A thorough description of the CTD-SRDL technology can be found in : 
  * Boehme L., Lovell P., Biuw M., Roquet F., Nicholson J., Thorpe S. E., Meredith M. P., and Fedak M., 2009. Technical Note: Animal-borne CTD-Satellite Relay Data Loggers for real-time oceanographic data collection. Ocean Sci., 5:685-695. doi: 10.5194/os-5-685-2009
* The delayed-mode general methodology and estimated accuracy of CTD-SRDL hydrographic data 
are presented in :
  * Roquet F., Charrassin J.-B., Marchand S., Boehme L., Fedak M., Reverdin G., and Guinet C., 2011. Validating hydrographic data obtained from seal-borne satellite-relayed data loggers. J. Atmos. Oceanic Technol., 28:787-801. doi: 10.1175/2010JTECHO801.1
* The density inversion removal algorithm is described in :
  * Barker, P. M. and McDougall, T. J., 2017. Stabilizing Hydrographic Profiles with Minimal Change to the Water Masses. J. Atmos. Oceanic Technol., 34:1935-1945. doi: 10.1175/JTECH-D-16-0111.1
* The thermal cell effect correction is described in :
  * Mensah, V., Roquet, F., Siegelman-Charbit, L., Picard, B., Pauthenet, E., Guinet, C., 2018.  A correction methodology for the thermal mass induced-errors of CTD tags mounted on marine mammals. J. Atmos. Oceanic Technol., 35:1237–1252. doi: 10.1175/JTECH-D-17-0141.1


### National specificities :
- For Australian data: 
Any users of IMOS data are required to clearly acknowledge the source of the material 
derived from IMOS in the format: 
"Data was sourced from the Integrated Marine Observing System (IMOS) - IMOS is a national 
collaborative research infrastructure, supported by Australian Government.” IMOS data is 
licensed under a Creative Commons Attribution (CCBY) License, 
(http://creativecommons.org.au/)."
- For German and South African data: 
Primary data are also made available through PANGAEA. Please cite :
doi10.1594/PANGAEA.150008 for data related to Marion Island (Southern Ocean Indian Sector)
doi10.1594/PANGAEA.150009 for data related to King George Island (Southern Ocean Atlantic Sector)
doi10.1594/PANGAEA.150010 for data related to Atka Bay, Drescher Inlet, Filchner Trough (Southern Ocean Atlantic Sector)    


