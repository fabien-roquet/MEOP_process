README: Important information about MEOP data products

Relevant to any version of the MEOP-CTD databases, MEOP-SMS databases and MEOP-TDR databases
Owner: Fabien Roquet (fabien.roquet@gu.se) and the MEOP consortium
Licence: https://opendatacommons.org/licenses/odbl/



BRIEF DESCRIPTION

Since 2004, several hundred thousands profiles of temperature and salinity have been 
collected by instrumented animals. The use of elephant seals has been particularly 
effective to sample the Southern Ocean and the North Pacific. These hydrographic data 
have been assembled in a quality-controlled database, the MEOP-CTD database, that can 
be accessed through this website.

The MEOP data portal distributes currently three different databases: 
  MEOP-CTD database: quality-controlled CTD profiles
  MEOP-SMS database: submesoscale-resolving high density CTD data
  MEOP-TDR database: high spatial density temperature/light data
  
For more information, visit the website meop.net
For any questions, contact info@meop.net or fabien.roquet@gu.se



DATA FORMATS

Binary netCDF formats for hydrographic profiles:
The ncARGO format: For a thorough scientific use of the data, or for oceanographic data centers, it is advised to use the marine mammal netCDF format (files in DATA_ncARGO) as it serves as the reference. This format can be easily read in Ocean Data View  using the Import/ARGO profiles/Float profiles menu, or using your favorite data processing software (e.g. Python, Matlab, IDL). Matlab tools and python tools are also available publicly to read and manipulate files in netCDF format. 

Binary netCDF formats for timeseries:
DATA_ncTRAJ: This is the format used to store and distribute timeseries data.

Text-file Formats :
DATA_csv_interp: A csv format (ASCII) is also provided (files in DATA_csv_interp) which can be opened with Excel or any text editor. Here, only data flagged as good are included, and are given on a regular vertical grid (1dbar spacing).
METADATA: Complementary text files compiling the names and values of metadata in a simple text file. The metadata found in this file consists of the list of global attributes and their values in the corresponding netCDF file.




HOW TO CITE

If you use this dataset for a publication, please add the following sentence in the Acknowledgement part:

"The marine mammal data were collected and made freely available by the International MEOP 
Consortium and the national programs that contribute to it (http://www.meop.net)."

Also please cite at least one of the following papers when you use the MEOP databases in your published work:
  Treasure, A.M., F. Roquet, I.J. Ansorge, M.N. Bester, L. Boehme, H. Bornemann, J.-B. Charrassin, D. Chevallier, D.P. Costa, M.A. Fedak, C. Guinet, M.O. Hammill, R.G. Harcourt, M.A. Hindell, K.M. Kovacs, M.-A. Lea, P. Lovell, A.D. Lowther, C. Lydersen, T. McIntyre, C.R. McMahon, M.M.C. Muelbert, K. Nicholls, B. Picard, G. Reverdin, A.W. Trites, G.D. Williams, and P.J.N. de Bruyn. 2017. Marine Mammals Exploring the Oceans Pole to Pole: A review of the MEOP consortium. Oceanography 30(2):132–138, doi: 10.5670/oceanog.2017.234
  Roquet F., Wunsch C., Forget G., Heimbach P., Guinet C., Reverdin G., Charrassin J.-B., Bailleul F., Costa D. P., Huckstadt L. A., Goetz K. T., Kovacs K. M., Lydersen C., Biuw M., Nøst O. A., Bornemann H., Ploetz, J., Bester M. N., Mcintyre T., Muelbert M. C., Hindell M. A., McMahon C. R., Williams G., Harcourt R., Field I. C., Chafik L., Nicholls K. W., Boehme L., and Fedak M. A., 2013. Estimates of the Southern Ocean General Circulation Improved by Animal-Borne Instruments. Geoph. Res. Letts., 40:1-5. doi: 10.1002/2013GL058304
  Roquet F., Williams G., Hindell M. A., Harcourt R., McMahon C. R., Guinet C., Charrassin J.-B., Reverdin G., Boehme L., Lovell P. and Fedak M. A., 2014. A Southern Indian Ocean database of hydrographic profiles obtained with instrumented elephant seals. Nature Scientific Data, 1:140028, doi: 10.1038/sdata.2014.28

Important technical papers about CTD-SRDL data: 
  Boehme L., Lovell P., Biuw M., Roquet F., Nicholson J., Thorpe S. E., Meredith M. P., and Fedak M., 2009. Technical Note: Animal-borne CTD-Satellite Relay Data Loggers for real-time oceanographic data collection. Ocean Sci., 5:685-695. doi: 10.5194/os-5-685-2009
  Roquet F., Charrassin J.-B., Marchand S., Boehme L., Fedak M., Reverdin G., and Guinet C., 2011. Validating hydrographic data obtained from seal-borne satellite-relayed data loggers. J. Atmosph. And Ocean. Tech., 28:787-801. doi: 10.1175/2010JTECHO801.1
  Siegelman, L., Roquet, F., Mensah, V., Rivière, P., Pauthenet, E., Picard, B., and Guinet, C., 2019. Correction and accuracy of high- and low-resolution CTD data from animal-borne instruments. Journal of Atmospheric and Oceanic Technology. doi: 10.1175/JTECH-D-18-0170.1

National specificities :
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



