# library(jsonlite)
library(tidyverse)

options(scipen = 999)


rm(list=ls())


source("Conversion Functions.R")
source("Auxiliary Functions.R")


# json_data <- fromJSON("jsonDefs.json")
# units_df <- json_data[[3]];rm(json_data)


# Loading the units from UCUN
load("units.Rdata")

units_df <- units_df %>% select(csCode_, ciCode_, magnitude_, dim_, cnv_,
                                cnvPfx_, isArbitrary_, moleExp_)

cnv_units <- units_df %>% filter(!is.na(cnv_))
arbitrary_units <- extractArbitrary(units_df)

# Specifying conversion units & values
amount <- 12
fromUnit <- "mg.%"
toUnit <- "kg.%"
molecularWeight <- 1.435



# Performing the conversion
convert(fromNum = amount, 
        fromUnit = fromUnit, 
        toUnit = toUnit,
        molecularWeight = molecularWeight)




# fromUnit <- "13.mg/dL3"
# toUnit <- "7.mmol2/L3"


# fromUnit <- "6/h.(1.g)/10^3.(11.h/(16.kg)).cm"
# toUnit <- "7/h.(12.g)/10^3.(3.h/(1.kg)).cm"


# fromUnit <- "Cel-3.10^3.16"
# toUnit <- "Cel-3.10^3.6"


# fromUnit <- "SIE/3"
# toUnit <- "mho.3"

# fromUnit <- "3.Cel.10^3"
# toUnit <- "8.Cel.10^3"



