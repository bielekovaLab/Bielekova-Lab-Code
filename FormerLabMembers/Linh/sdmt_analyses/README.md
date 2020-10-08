# Hello! Thank you for reading our smartphone SDMT paper!
And for coming by to check out our data and scripts. Below are some information to help you get started.

**If you are a reviewer, please email Vanessa Morgan at vanessa.morgan@nih.gov for access to the raw data files. The raw files will be added for public view after the paper is published. Thank you**

You will need R to run all of these scripts. The best way to get everything is to clone the repository. Follow this tutorial if you've never cloned a repository before (https://www.youtube.com/watch?v=YxZ8J2rqhEM).

This was run on R version 3.5.1 (2018-07-02). Package versions can be found in the paper's reference 12-31.

IMPORTANT: For training and validation cohorts division, set.seed() works differently for R version 3.6.0 and above (https://stackoverflow.com/questions/47199415/is-set-seed-consistent-over-different-versions-of-r-and-ubuntu). If you are running an R version that is at least 3.6.0, be sure to run the following line in your console before doing anything else to get reproducible results:

RNGkind(sample.kind = "Rounding")

# Data
Raw data (not cleaned for outliers and duplicates, etc.) for the SDMT tests, neurological disability scores, hand tapping data, and MRI data are located in the "data_raw" folder. Data in the "data_raw" folder serve as input for the data cleaning scripts, located in the "cleaning_scripts" folder. Outputs of the cleaning scripts go out to the "data_clean" folder. Everything in the "data_clean" folder then get used in the analyses scripts.
Any SDMT data in the "data_raw" folder do not have accurate age. Use the raw data in "data_raw_age" to get accurate age. 

To navigate the "cleaning_scripts" folder, first go to the script "AllClean_SDMT" that is in that folder. This script utilizes the scripts "SymbolCleanFuncPreV1" and "SymbolCleanFuncV1", as well as all of the scripts that are located in the "sdmt_clean" folder. Running the AllClean script will clean ALL of the SDMT data (app and written).

Data cleaning for everything else are labeled in the format "data_clean" (i.e. mri_clean, neurex_clean, etc.)

I tried my best to put comment in the data cleaning scripts so that you can understand both the scripts and the structure of the raw input data. Please reach out to the corresponding author if you have any questions.

# Analyses scripts
All cross-sectional analyses that are documented in the paper can be found in the script "cross_analyses"; longitudinal analyses can be found in "long_analyses"; analyses relating to redefining the clinically significant change threshold for the written SDMT (figure 10 in the paper) are located in "writ_change".

The functions needed to run these scripts are in the folders "funcs_cross" and "funcs_long". The "writ_change" script requires one function in the "funcs_cross" folder and the "app_match_long" function. See comments in all of the analyses scripts for more details.

Analyses relating to supplemental figures for cross-sectional and longitudinal analyses are at the bottom of the cross and long analyses scripts.

#Results
Results from cross-sectional and longitudinal analyses can be viewed directly in R (from running the analyses scripts), but if you would like to see the output images, you can find them in the "results_cross" and "results_long" folders. Results from the "writ_change" script will also be in the "results_long" folder.

Demographics information (from table 1 and 2 in the paper) can be grabbed using the demo_grab script.

# Contact info
Reach out to Dr. Bibi Bielekova (bibi.bielekova@nih.gov) for any questions and most up-to-date contact information for Linh (repository's author).
