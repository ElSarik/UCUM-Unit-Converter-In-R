
#===============================================#
####          Conversion Function            ####
#===============================================#

convert <- function(fromNum, fromUnit, toUnit, molecularWeight){
  
  # Checking if fromNum is not numeric
  if(!is.numeric(fromNum)){
    stop(paste0("Amount must be numeric! Given value: ", fromNum))
  }
  
  # Checking if fromUnit exists in units_df
  from_unit <- check_unit_existence(fromUnit, 
                                    parsed_unit = FALSE)
  
  # If the result is NULL, parse the unit and add it to units_df
  if(is.null(from_unit)){
    new_unit <- parseUnit(fromUnit)
    
    print(new_unit)
    
    units_df <<- rbind(units_df, new_unit)
    
    # Checking if fromUnit exists in units_df
    # This time it should always exist (if valid) 
    # as it has just been added to the units_df
    from_unit <- check_unit_existence(fromUnit, 
                                      parsed_unit = FALSE)
  }
  
  fromCsCode_ <- from_unit$csCode_
  fromIsArbitrary_ <- from_unit$isArbitrary_
  fromDimension_ <- eval(parse(text = from_unit$dim_$dimVec_))
  
  
  # Checking if toUnit exists in units_df
  to_unit <- check_unit_existence(toUnit, 
                                  parsed_unit = FALSE)
  
  # If the result is NULL, parse the unit and add it to units_df
  if(is.null(to_unit)){
    new_unit <- parseUnit(toUnit)
    
    print(new_unit)
    
    units_df <<- rbind(units_df, new_unit)
    
    # Checking if toUnit exists in units_df
    # This time it should always exist (if valid) 
    # as it has just been added to the units_df
    to_unit <- check_unit_existence(toUnit, 
                                      parsed_unit = FALSE)
  }
  
  toCsCode_ <- to_unit$csCode_
  toIsArbitrary_ <- to_unit$isArbitrary_
  toDimension_ <- eval(parse(text = to_unit$dim_$dimVec_))
  
  
  # print(fromCsCode_)
  # print(fromIsArbitrary_)
  # print(fromDimension_)
  # 
  # print(toCsCode_)
  # print(toIsArbitrary_)
  # print(toDimension_)
  
  
  
  # Checking if we are converting from/to arbitrary units
  if(isTRUE(toIsArbitrary_)){
    stop(paste0("Attempt to convert to arbitrary unit ", toCsCode_))
  }
  if(isTRUE(fromIsArbitrary_)){
    stop(paste0("Attempt to convert arbitrary unit ", fromCsCode_))
  }
  
  # Checking whether the given units are Mole-To-Mass commensurable
  if(!isTRUE(all.equal(fromDimension_, toDimension_))){
    isMMcom <- isMoleMassCommensurable(from_unit, to_unit)
    
    isMMcom_status <- isMMcom$status
    isMMcom_function <- isMMcom$fun
    
    if(isTRUE(isMMcom_status) & !is.numeric(molecularWeight)){
      stop("A numeric molecular weight must be provided!")
    }
    
    if(!isTRUE(isMMcom_status)){
      stop(paste0("Sorry. ", fromCsCode_, " cannot be converted to ", toCsCode_))
    }
    
    if(isTRUE(isMMcom_status) & is.numeric(molecularWeight)){
      if(isMMcom_function == 'MassToMole'){
        result <- convertMassToMol(fromNum, fromUnit, toUnit, molecularWeight)
        return(result)
      }
      if(isMMcom_function == 'MoleToMass'){
        result <- convertMolToMass(fromNum, fromUnit, toUnit, molecularWeight)
        return(result)
      }
    }
  }
  
  result <- convertFrom(fromNum, fromUnit, toUnit)
  return(result)
}




#===============================================#
####          Regular Conversions            ####
#===============================================#


# Performs a regular conversion between 2 units
convertFrom <- function(fromNum, fromUnit, toUnit){
  
  # filtering the units and extracting the necessary parameters of the
  # 'from' unit
  from_unit <- check_unit_existence(fromUnit, 
                                    parsed_unit = TRUE)
  
  fromCnv <- from_unit$cnv_
  fromMag <- from_unit$magnitude_
  fromCnvPfx <- from_unit$cnvPfx_
  
  # filtering the units and extracting the necessary parameters of the
  # 'to' unit
  to_unit <- check_unit_existence(toUnit, 
                                  parsed_unit = TRUE)
  
  toCnv <- to_unit$cnv_
  toMag <- to_unit$magnitude_
  toCnvPfx <- to_unit$cnvPfx_
  
  
  # Checking if 'fromCnv' has a value or not ('cnv_' column in units data)
  # If it contains a value, we perform the specific conversion of the value
  # located in 'convert_from_to()' function.
  if(!is.na(fromCnv)){
    
    # We get the "from" conversion of the given 'fromCnv' value
    from_x_conversion <- convert_from_to(fromCnv = fromCnv,
                                         x = fromNum * fromCnvPfx,
                                         conv_status = "FROM")
    
    # We multiply the received value with the magnitude of 'fromUnit'
    from_x <- from_x_conversion * fromMag
  }
  
  # If no specific conversion is required, we just multiply the given
  # value (fromNum) with the magnitude of 'fromUnit'
  else{
    from_x <-fromNum * fromMag
  }
  
  
  # At this point the given value (fromNum) has now either been converted
  # to a base form (as specified in 'convert_from_to()' function), or it has
  # been multiplied with the magnitude of 'fromUnit'
  
  # The same procedure is now being performed with the new value, but for
  # the 'toUnit' conversion
  if(!is.na(toCnv)){
    to_x_conversion <- convert_from_to(fromCnv = toCnv,
                                       x = from_x / toMag,
                                       conv_status = "TO")
    
    to_x <- to_x_conversion / toCnvPfx
  }
  else{
    to_x <- from_x / toMag
  }
  
  return(to_x)
}


#===============================================#
####          Convert Mass to Mole           ####
#===============================================#

# Converts from a Mass unit to a Mole unit
convertMassToMol <- function(fromNum, fromUnit, toUnit, molecularWeight){
  # 'fromUnit' corresponds to the Mass unit
  # 'toUnit' corresponds to the Mole unit
  
  # Extracting magnitude of the 'from' unit
  from_unit <- check_unit_existence(fromUnit, 
                                    parsed_unit = TRUE)
  
  fromMag <- from_unit$magnitude_
  
  # Extracting magnitude of the 'to' unit
  to_unit <- check_unit_existence(toUnit, 
                                  parsed_unit = TRUE)
  
  toMag <- to_unit$magnitude_
  
  # Extracting magnitude of base 'mole' unit
  mol_unit <- units_df %>% filter(csCode_ == 'mol')
  
  molMag <- mol_unit$magnitude_
  
  
  # The prefix values that have been applied to this unit, which is the mass
  # (grams) unit, are reflected in the magnitude.  So the number of moles
  # represented by this unit equals the number of grams -- amount * magnitude
  # divided by the molecular Weight
  molAmount <- (fromMag * fromNum) / molecularWeight
  
  # We acquire the base mole magnitude and divide it out of the current magnitude.
  molesFactor <- toMag / molMag
  
  # return the molAmt divided by the molesFactor as the number of moles
  # for the molUnit
  return(molAmount / molesFactor)
}



#===============================================#
####          Convert Mole to Mass           ####
#===============================================#

# Converts from a Mole unit to a Mass unit
convertMolToMass <- function(fromNum, fromUnit, toUnit, molecularWeight){
  # 'fromUnit' corresponds to the Mole unit
  # 'toUnit' corresponds to the Mass unit
  
  # Extracting magnitude of the 'from' unit
  from_unit <- check_unit_existence(fromUnit, 
                                    parsed_unit = TRUE)
  
  fromMag <- from_unit$magnitude_
  
  # Extracting magnitude of the 'to' unit
  to_unit <- check_unit_existence(toUnit, 
                                  parsed_unit = TRUE)
  
  toMag <- to_unit$magnitude_
  
  # Extracting magnitude of base 'mole' unit
  mol_unit <- units_df %>% filter(csCode_ == 'mol')
  
  molMag <- mol_unit$magnitude_
  
  
  # Determine what prefix values (mg or mg/dL, etc.) have been applied to
  # this unit by dividing the simple mole unit magnitude out of the
  # current mole unit magnitude.
  molesFactor <- fromMag / molMag
  
  # The number of grams (mass) is equal to the number of moles (fromNum)
  # times the molecular weight.  We also multiply that by the prefix values
  # applied to the current unit (molesFactor) to get the grams for this
  # particular unit.
  massAmount <- (molesFactor * fromNum) * molecularWeight
  
  # Finally, we return the mass amount/grams for this particular unit
  # divided by any effects of prefixes applied to the "to" unit, which
  # is assumed to be some form of a gram unit
  return(massAmount / toMag)
}



#===============================================#
####          Specific Conversions           ####
#===============================================#

# Based on the conversions from here:
# https://github.com/lhncbc/ucum-lhc/blob/a8d3182db3a83b0a884c5a4988c44ca7e1330b83/source/ucumFunctions.js

# All entries of the "cnv_" column of "jsonDefs.json" should be included!

convert_from_to <- function(fromCnv, x, conv_status = 'FROM'){

  # // Celsius - convert to Celsius from kelvin and from Celsius to kelvin
  # // where kelvin is the base unit for temperature
  if(fromCnv == "Cel"){
    if(conv_status == 'TO'){
      return(x - 273.15)
    }
    if(conv_status == 'FROM'){
      return(x + 273.15)
    }
  }
  
  
  # // Fahrenheit - convert to Fahrenheit from kelvin and from Fahrenheit to
  # // kelvin - which is the base unit for temperature
  if(fromCnv == "degF"){
    if(conv_status == 'TO'){
      return(x - 459.67)
    }
    if(conv_status == "FROM"){
      return(x + 459.67)
    }
  }
  
  
  # // Reaumur - convert between Reaumur and Kelvin.   Because of the way the
  # // calling code in the Units class is set up (in the convertFrom method),
  # // what is given here as the convertTo function is actually the convert
  # // from method and vice versa
  if(fromCnv == "degRe"){
    if(conv_status == 'TO'){
      return(x - 273.15)
    }
    if(conv_status == "FROM"){
      return(x + 273.15)
    }
  }
  
  
  # // pH - convert to pH from moles per liter and from moles per liter to pH
  # // where a mole is an amount of a substance (a count of particles)
  if(fromCnv == "pH"){
    if(conv_status == 'TO'){
      return(- log(x) / log(10))
    }
    if(conv_status == "FROM"){
      return(10^-x)
    }
  }
  
  
  # // ln - natural logarithm (base e 2.71828) - apply (cnvTo) and invert (cnvFrom)
  # // and 2ln - two times the natural logarithm
  if(fromCnv == "ln"){
    if(conv_status == 'TO'){
      return(log(x))
    }
    if(conv_status == "FROM"){
      return(exp(x))
    }
  }
  
  if(fromCnv == "2ln"){
    if(conv_status == 'TO'){
      return(2 * log(x))
    }
    if(conv_status == "FROM"){
      return(exp(x / 2))
    }
  }
  
  
  # // lg - the decadic logarithm (base 10)
  if(fromCnv == "lg"){
    if(conv_status == 'TO'){
      return(log(x) / log(10))
    }
    if(conv_status == "FROM"){
      return(10^x)
    }
  }
  
  if(fromCnv == "10lg"){
    if(conv_status == 'TO'){
      return(10 * log(x) / log(10))
    }
    if(conv_status == "FROM"){
      return(10^(x/10))
    }
  }
  
  if(fromCnv == "20lg"){
    if(conv_status == 'TO'){
      return(20 * log(x) / log(10))
    }
    if(conv_status == "FROM"){
      return(10^(x/20))
    }
  }
  
  # // The plain text ucum units file uses 'lgTimes2'
  if(fromCnv == "lgTimes2"){
    if(conv_status == 'TO'){
      return(2 * log(x) / log(10))
    }
    if(conv_status == "FROM"){
      return(10^(x/2))
    }
  }
  
  
  # // ld - dual logarithm (base 2)
  if(fromCnv == "ld"){
    if(conv_status == 'TO'){
      return(log(x) / log(2))
    }
    if(conv_status == "FROM"){
      return(2^x)
    }
  }
  
  
  # // tan - tangent
  if(fromCnv == "tanTimes100"){
    if(conv_status == 'TO'){
      return(tan(x) * 100)
    }
    if(conv_status == "FROM"){
      return(atan(x / 100))
    }
  }
  if(fromCnv == "100tan"){
    if(conv_status == 'TO'){
      return(tan(x) * 100)
    }
    if(conv_status == "FROM"){
      return(atan(x / 100))
    }
  }
  
  
  
  # // sqrt - square root
  if(fromCnv == "sqrt"){
    if(conv_status == 'TO'){
      return(sqrt(x))
    }
    if(conv_status == "FROM"){
      return(atan(x * x))
    }
  }
  
  
  
  # // inv - inverse
  if(fromCnv == "inv"){
    if(conv_status == 'TO'){
      return(1.0 / x)
    }
    if(conv_status == "FROM"){
      return(1.0 / x)
    }
  }
  
  
  # // homeopathic potency functions
  if(fromCnv == "hpX"){
    if(conv_status == 'TO'){
      return(- (log(x) / log(10)))
    }
    if(conv_status == "FROM"){
      return(10^-x)
    }
  }
  
  if(fromCnv == "hpC"){
    if(conv_status == 'TO'){
      return(- (log(x) / log(100)))
    }
    if(conv_status == "FROM"){
      return(100^-x)
    }
  }
  
  if(fromCnv == "hpM"){
    if(conv_status == 'TO'){
      return(- (log(x) / log(1000)))
    }
    if(conv_status == "FROM"){
      return(1000^-x)
    }
  }
  
  
  if(fromCnv == "hpQ"){
    if(conv_status == 'TO'){
      return(- (log(x) / log(50000)))
    }
    if(conv_status == "FROM"){
      return(50000^-x)
    }
  }
}


