# PCB water concentration project

## License

PCBWaterProject is licensed under the 2-Clause BSD License - see the [LICENSE](LICENSE) file for details.

----------------------
General Information
----------------------

Deposit Title: Water concentration of PCBs from surfaces water from USA

Contributor information:

Andres Martinez, PhD
University of Iowa - Department of Civil & Environmental Engineering
Iowa Superfund Research Program (ISRP)
andres-martinez@uiowa.edu
ORCID: 0000-0002-0572-1494

This README file was generated on September 24, 2024 by Andres Martinez.

This work was supported by the National Institutes of Environmental Health Sciences (NIEHS) grant #P42ES013661.  The funding sponsor did not have any role in study design; in collection, analysis, and/or interpretation of data; in creation of the dataset; and/or in the decision to submit this data for publication or deposit it in a repository.

This README file describes the codes generated to analyse and map water PCB concentrations from Superfund sites used in this paper:

The script is developed for the analysis and mapping of water concentrations of PCBs from Superfund sites,Project 4, Aim 3, ISRP.

--------
PREREQUISITES & DEPENDENCIES
--------

This section of the ReadMe file lists the necessary software required to run codes in "R".

Software:
- Any web browser (e.g., Google Chrome, Microsoft Edge, Mozilla Firefox, etc.)
- R-studio for easily viewing, editing, and executing "R" code as a regular "R script" file:
https://www.rstudio.com/products/rstudio/download/

--------
SOFTWARE INSTALLATION
--------

This section of the ReadMe file provides short instructions on how to download and install "R Studio".  "R Studio" is an open source (no product license required) integrated development environment (IDE) for "R" and completely free to use.  To install "R Studio" follow the instructions below:

1. Visit the following web address: https://www.rstudio.com/products/rstudio/download/
2. Click the "download" button beneath RStudio Desktop
3. Click the button beneath "Download RStudio Desktop".  This will download the correct installation file based on the operating system detected.
4. Run the installation file and follow on-screen instructions. 

--------
R FILES AND ESTRUCTURE
--------
It is recommended to create a Project in R (e.g., PCBWaterProject.Rproj). Download the project file and the R folder where the codes are, and the Subfolder.R code. Run first the Subfolder code, which will generate all the subfolders. 

--------
DATA
--------

PCB data were compiled from various sources and processed both manually and using R. Examples of how the data were extracted and wrangled can be found in the 'R/ExtractingData/' directory.

The final dataset is available at:
Martinez, Andres (2024): Dataset of surface water concentrations of Polychlorinated Biphenyls in the U.S. from 1979–2020 [dataset]. PANGAEA, https://doi.org/10.1594/PANGAEA.972705

The dataset can also be downloaded using the script found at 'R/Pangaea/PangaeaDownloadDataset.R'.

Water temperatures and flows were mostly obtained from USGS station, uisng the package in R dataRetrieval and library dataRetrieval.





