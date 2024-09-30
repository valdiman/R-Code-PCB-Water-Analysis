# R code to download and format data from Pangaea

# Packages and libraries needed -------------------------------------------------------------------
# Install packages
install.packages('pangaear')
# Library
library(pangaear) # Read data from Pangaea

# Read data from Pangaea and format data ------------------------------

# Read data from Pangaea repository
# Citation:
# Martinez, Andres (publication year): Dataset of surface water
# concentrations of Polychlorinated Biphenyls in the U.S.
# from 1979 - 2020 [dataset].  PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.972705
# final doi: https://doi.org/10.1594/PANGAEA.972705

# Set cache path to the project folder
pg_cache$cache_path_set(full_path = "Data/")
# Download original datasets from Pangaea
d.0 <- pg_data(doi = '10.1594/PANGAEA.972705')

# Extract dataset
d <- d.0[[1]]$data

# Need to change the columns name to they match the R codes
# Vector with all the new column names provided
new_col_names <- c(
  "Source", "SampleID", "EPARegion", "StateSampled", "LocationName", 
  "SiteName", "SiteID", "SampleDate", "Latitude", "Longitude", 
  "PhaseMeasured", "EPAMethod", "AroclorCongener", "PCB1", "PCB2", "PCB3", 
  "PCB4+10", "PCB5+8", "PCB6", "PCB7+9", "PCB11", "PCB12+13", "PCB14", 
  "PCB15", "PCB16+32", "PCB17", "PCB18+30", "PCB19", 
  "PCB20+21+28+31+33+50+53", "PCB22", "PCB23", "PCB24+27", "PCB25", 
  "PCB26+29", "PCB34", "PCB35", "PCB36", "PCB37+42", "PCB38", "PCB39", 
  "PCB40+41+64+71+72", "PCB43+49+52+69+73", "PCB44+47+65", "PCB45+51", 
  "PCB46", "PCB48+59+62+75", "PCB54", "PCB55", "PCB56+60", "PCB57", 
  "PCB58", "PCB61+66+70+74+76+93+95+98+100+102", "PCB63", "PCB67", 
  "PCB68", "PCB77+85+110+111+115+116+117", "PCB78", "PCB79", "PCB80", 
  "PCB81+86+87+97+107+108+109+112+119+124+125", 
  "PCB82+135+144+151+154", "PCB83+99", "PCB84+92", "PCB88+91", "PCB89", 
  "PCB90+101+113", "PCB94", "PCB96", "PCB103", "PCB104", "PCB105", 
  "PCB106+118", "PCB114+122+131+133+142+146+165", "PCB120", "PCB121", 
  "PCB123+139+140+147+149", "PCB126", "PCB127", "PCB128+162+166+167", 
  "PCB129+137+138+158+160+163+164+176+178", "PCB130", 
  "PCB132+153+161+168", "PCB134+143", "PCB136", "PCB141", "PCB145", 
  "PCB148", "PCB150", "PCB152", "PCB155", "PCB156+157+172+197+200", 
  "PCB159", "PCB169", "PCB170+190", "PCB171+173+202", "PCB174", 
  "PCB175", "PCB177", "PCB179", "PCB180+193", "PCB181", 
  "PCB182+187", "PCB183+185", "PCB184", "PCB186", "PCB188", "PCB189", 
  "PCB191", "PCB192", "PCB194+205", "PCB195+208", "PCB196+203", 
  "PCB198+199+201", "PCB204", "PCB206", "PCB207", "PCB209", "A1016", 
  "A1221", "A1232", "A1242", "A1248", "A1254", "A1260", "tPCB"
)

# Apply the new column names
colnames(d) <- new_col_names

# Export modified dataset
write.csv(d, file = "Data/USAWaterPCB.csv",
          row.names = FALSE)





