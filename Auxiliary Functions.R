# Checking whether the given unit exists inside the units_df table
# and returning the found unit.

# parsed_unit indicates that the given unit is parsed, and not a
# unit expression.
check_unit_existence <- function(unit, parsed_unit = TRUE){
  
  # Checking if unit exists as csCode_
  found_unit <- units_df %>% filter(csCode_ == unit)
  
  # No match based on given csCode_
  if(nrow(found_unit) < 1){
    
    # Checking if unit exists as ciCode_
    found_unit <- units_df %>% filter(ciCode_ == toupper(unit))
    
    # No match based on given ciCode_
    if(nrow(found_unit) < 1){
      
      if(isTRUE(parsed_unit)){
        stop(paste("The unit does not exist: ", unit))
      }
      else{
        return(NULL)
      }
    }
    # return the found unit based on the ciCode_
    else{
      return(found_unit)
    }
  }
  # return the found unit based on the csCode_
  else{
    return(found_unit)
  }
}






# Extracting arbitrary units in a format suitable to detect
# whether a unit expression is arbitrary or not
extractArbitrary <- function(units_df){
  
  arbitrary_units <- units_df %>% filter(isArbitrary_ == TRUE)
  
  arbitrary_units <- arbitrary_units$csCode_
  
  arb_units_fixed <- c()
  for(unit in arbitrary_units){
    format_unit <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", unit, perl = TRUE)
    
    arb_units_fixed <- c(arb_units_fixed, format_unit)
  }
  
  arb_units_fixed <- paste(arb_units_fixed, collapse = "|")
  return(arb_units_fixed)
}



# Attempting to parse the unit
parseUnit <- function(uStr){
  
  # Removing any white spaces inside the string
  uStr <- gsub(" ", "", uStr)
  
  # Replacing "{...}" with 1 or ""
  if(str_detect(uStr, "\\w+\\{[^{}]+\\}")){
    uStr <- str_replace(uStr, "\\{[^{}]+\\}", "")
  }
  uStr <- str_replace(uStr, "\\{[^{}]+\\}", "1")
  
  # Checking if a unit ends in an operator ("." or "/")
  if(str_ends(uStr, "\\.") | str_ends(uStr, "\\/")){
    stop(paste0(uStr, " is not a valid unit expression! \n
                Unit expressions cannot end with an operator"))
  }
  
  # Checking if a unit contains multiple operators in a row ("." or "/")
  if(grepl("[\\/\\.]{2}|[\\/\\.][\\/\\.]", uStr)){
    stop(paste0(uStr, " is not a valid unit expression! \n
                 Unit expressions cannot contain multiple operators in a row"))
  }
  
  # Checking if all parentheses inside the string are valid
  if(!isTRUE(check_parentheses(uStr))){
    stop(paste0(uStr, " is not a valid unit expression! \n
                Unit expression contains invalid parentheses"))
  }
  
  # Checking if addition (+) or subtraction (-) operations are being performed
  if(str_detect(uStr, "\\d+[+|-]+\\d+")){
    stop(paste0(uStr, " is not a valid unit expression! \n
                 Unit expressions cannot contain +/- operations"))
  }
  
  
  #!!!!!!!!!!!!!!!!! CHECK THIS AGAIN LATER !!!!!!!!!!!!!!!!#
  # ========================================================#
  # # Checking if parentheses close and open without a valid expression in between
  # if(!grepl("(\\)[.|\\/].+?[.|\\/]\\()", uStr) & !grepl("(\\)[.|\\/]\\()", uStr)){
  #   print(paste0(uStr, " is not a valid unit expression!"))
  #   print("A valid unit expression must be provided between ')' and '(' ")
  # }
  
  
  original_uStr <- uStr
  
  # Fixing exponents and multiplication symbols
  uStr <- str_replace_all(uStr, '\\*', '^')
  uStr <- str_replace_all(uStr, '\\.', '*')
  uStr <- processExponents(uStr)
  
  # Adding 1 in front of unit expression if it starts with "/"
  if(str_starts(uStr, "\\/")){
    uStr <- paste0("1", uStr)
  }
  
  # Checking if unit expression contains an exponent other than
  # valid unit exponents (eg h^3, 10^4)
  if(str_detect(uStr, "\\^")){
    if(isTRUE(str_detect(uStr, "\\b(?!10\\^)(\\d+\\^\\d+)\\b"))){
      stop(paste0(original_uStr, " is not a valid unit expression! \n
                   Unit expressions may only contain valid unit exponents"))
    }
  }
  
  
  # Checking if uStr contains a non-ratio unit, which require
  # special treatment
  cnv_result <- identify_cnv(uStr)
  
  cnv <- cnv_result$cnv_exists
  cnv_csCode_ <- cnv_result$cs_code
  cnv_value <- cnv_result$cnv_value
  
  
  # Checking if unit expression contains an arbitrary unit
  if(str_detect(uStr, arbitrary_units)){
    arbitrary <- TRUE
  }
  else{
    arbitrary <- FALSE
  }
  
  
  if(isTRUE(cnv)){
    
    # Calculate the dimension of the given unit expression
    unit_dimension <- calculate_dimension(uStr)
    
    unit_dimension <- paste0("c(",paste(unit_dimension, collapse = ", "),")")
    
    # Calculate the mole exponent of the given unit expression
    moleExp_ <- calculate_moleExponent(uStr)
    
    # Calculating cnv_ magnitude and Pfx_
    cnv_mag_pfx <- calculate_cnv_mag_pfx(uStr, cnv_csCode_)
    unit_mag <- cnv_mag_pfx$cnv_magnitude
    cnvPfx_ <- cnv_mag_pfx$cnvPfx_
    
    # print(uStr)
    # print(unit_mag)
    # print(unit_dimension)
    # print(moleExp_)
    # print(cnv_value)
    # # print(cnv_csCode_)
    # print(cnvPfx_)
    # print(arbitrary)
    
    new_unit <- data.frame(csCode_ = original_uStr,
                           ciCode_ = original_uStr,
                           magnitude_ = unit_mag,
                           dim_ = unit_dimension,
                           cnv_ = cnv_value,
                           cnvPfx_ = cnvPfx_,
                           isArbitrary_ = arbitrary,
                           moleExp_ = moleExp_)
    
    
    return(new_unit)
    
  }else{
    
    # Processing the unit expression by separating the parentheses
    uStr_processed <- processParens(uStr, check_nested_parentheses(uStr))
    # print(uStr_processed)
    
    # Processing the unit expression by separating the operators from the
    # beginning and end of each string element
    uStr_processed <- processOperators(uStr_processed)
    # print(uStr_processed)
    
    # Calculating the magnitude of the given unit expression
    unit_mag <- calculate_magnitude(uStr_processed)
    
    # Calculate the dimension of the given unit expression
    unit_dimension <- calculate_dimension(uStr)
    
    unit_dimension <- paste0("c(",paste(unit_dimension, collapse = ", "),")")
    
    # Calculate the mole exponent of the given unit expression
    moleExp_ <- calculate_moleExponent(uStr)
    
    # cnvPfx_ is set to 1
    cnvPfx_ <- 1
    
    # print(uStr)
    # print(unit_mag)
    # print(unit_dimension)
    # print(moleExp_)
    # print(cnv_value)      # cnv_value should be NA
    # print(cnvPfx_)
    # print(arbitrary)
    
    
    new_unit <- data.frame(csCode_ = original_uStr,
                           ciCode_ = original_uStr,
                           magnitude_ = unit_mag,
                           dim_ = unit_dimension,
                           cnv_ = cnv_value,
                           cnvPfx_ = cnvPfx_,
                           isArbitrary_ = arbitrary,
                           moleExp_ = moleExp_)
    
    
    return(new_unit)
  }
}




# Checking if the 2 given units are mole-to-mass commensurable
isMoleMassCommensurable <- function(fromUnit, toUnit){
  
  fromDimension_ <- eval(parse(text = fromUnit$dim_$dimVec_))
  fromMoleExp_ <- fromUnit$moleExp_
  
  toDimension_ <- eval(parse(text = toUnit$dim_$dimVec_))
  toMoleExp_ <- toUnit$moleExp_
  
  
  if(toMoleExp_ == 1 & fromMoleExp_ == 0){
    
    test_dim <- toDimension_
    
    test_dim[3] <- test_dim[[3]] + toMoleExp_
    
    
    if(isTRUE(all.equal(test_dim, fromDimension_))){
      return(list('status' = TRUE,
                  'fun' = 'MassToMole'))
    }
    else{
      return(list('status' = FALSE,
                  'fun' = 'MassToMole'))
    }
  }
  
  if(fromMoleExp_ == 1 & toMoleExp_ == 0){
    
    test_dim <- fromDimension_
    
    test_dim[3] <- test_dim[3] + fromMoleExp_
    
    if(isTRUE(all.equal(test_dim, toDimension_))){
      return(list('status' = TRUE,
                  'fun' = 'MoleToMass'))
    }
    else{
      return(list('status' = FALSE,
                  'fun' = 'MoleToMass'))
    }
  }
  
  return(list('status' = FALSE,
              'fun' = NA))
}



# Checking if the parentheses of a string (if any) are balanced:
# - All parentheses open and close
# - No parentheses close before being opened
# - No parentheses remain opened
check_parentheses <- function(input_string) {
  
  # Initialize a counter to keep track of open parentheses
  open_parentheses_count <- 0
  
  # Iterate through each character in the input string
  for (char in strsplit(input_string, "")[[1]]) {
    if (char == "(") {
      # Increment the counter if an open parenthesis is encountered
      open_parentheses_count <- open_parentheses_count + 1
    } else if (char == ")") {
      # Decrement the counter if a close parenthesis is encountered
      open_parentheses_count <- open_parentheses_count - 1
      
      # If the counter goes negative, there was a closing parenthesis without an opening one
      if (open_parentheses_count < 0) {
        return(FALSE)
      }
    }
  }
  
  # If the counter is not zero at the end, there were more open than close parentheses
  return(open_parentheses_count == 0)
}



# Checking if a unit string contains nested parentheses or not,
# as well as whether it contains unclosed parentheses
check_nested_parentheses <- function(input_string) {
  stack <- 0
  max_depth <- 0
  
  for (char in strsplit(input_string, "")[[1]]) {
    if (char == "(") {
      stack <- stack + 1
      if (stack > max_depth) {
        max_depth <- stack
      }
    } else if (char == ")") {
      if (stack == 0) {
        return(FALSE)  # Closing parenthesis without a corresponding opening parenthesis
      } else {
        stack <- stack - 1
      }
    }
  }
  
  if (stack > 0) {
    return(FALSE)  # Unclosed parentheses
  } else if (max_depth > 1) {
    return(TRUE)   # Nested parentheses
  } else {
    return(FALSE)  # No nested parentheses
  }
}



# Separating the input string based on the parentheses it contains.
# A list is returned with each parenthesis segment being in a separate entry
# inside the list.
# eg. uStr <- "g/(11.h/(16.kg)).cm" gets converted to:
#
# "g/"
# "(11.h/"
# "(16.kg))"
# ".cm"
processParens <- function(uStr, nested_parens = FALSE){
  result <- c()
  character_pos <- 0
  character_start <- 1
  
  open_paren_pos <- 0
  open_paren_status <- 0
  
  close_paren_pos <- 0
  close_paren_status <- 0
  
  # The string contains nested parentheses
  if(isTRUE(nested_parens)){
    for(character in strsplit(uStr, "")[[1]]){
      
      # print(character)
      
      character_pos <- character_pos + 1
      
      # If an open parenthesis is found, add all characters until the parenthesis
      # to the result list.
      if(character == "("){
        open_paren_pos <- character_pos
        
        str_before_paren <- substr(uStr, character_start, open_paren_pos - 1)
        result <- c(result, str_before_paren)
        
        character_start <- character_pos
      }
      
      # If a closed parenthesis is found, check if the following character is
      # also a closed parenthesis. If yes, do nothing. if no, add all characters
      # until the parenthesis to the result list
      if(character == ")"){
        
        if(length(strsplit(uStr, "")[[1]]) != character_pos){
          if(strsplit(uStr, "")[[1]][character_pos + 1] != ")"){
            close_paren_pos <- character_pos
            
            str_before_paren <- substr(uStr, character_start, close_paren_pos)
            result <- c(result, str_before_paren)
            
            character_start <- character_pos + 1
          }
        }
      }
      
      # If character_pos has reached the end of the string
      # add all the remaining characters and add it to the result list.
      if(character_pos == length(strsplit(uStr, "")[[1]])){
        str_remaining <- substr(uStr, character_start, character_pos)
        result <- c(result, str_remaining)
      }
    }
  }
  
  # The string contains regular parentheses
  else{
    for(character in strsplit(uStr, "")[[1]]){
      character_pos <- character_pos + 1
      
      if(character == "("){
        open_paren_pos <- character_pos
        open_paren_status <- 1
      }
      
      if(character == ")"){
        close_paren_pos <- character_pos
        close_paren_status <- 1
      }
      
      # Once a parenthesis has opened and closed, add all characters before
      # the parenthesis to the results list
      # Then add the parenthesis to the results list, reset and continue
      # looking through the characters.
      if(open_paren_status == 1  & close_paren_status == 1){
        str_before_paren <- substr(uStr, character_start, open_paren_pos - 1)
        parenthesis <- substr(uStr, open_paren_pos, close_paren_pos)
        
        result <- c(result, str_before_paren, parenthesis)
        character_start <- close_paren_pos + 1
        open_paren_status <- 0
        close_paren_status <- 0
      }
    }
    
    if(close_paren_pos != character_pos){
      str_after_last_paren <- substr(uStr, close_paren_pos + 1, character_pos)
      result <- c(result, str_after_last_paren)
    }
  }
  
  return(as.list(result))
}


# Splitting the elements of the given list based on whether an operator
# exists either as the first or last character of the element string
processOperators <- function(uStr_list){
  
  result <- c()
  
  for(element in uStr_list){
    
    # Making sure that the string element is larger than 1
    if(length(strsplit(element, "")[[1]]) > 1){
      
      # Checking if string element starts AND ends with an operator
      if((str_starts(element, "\\*") | str_starts(element, "\\/")) &
         (str_ends(element, "\\*") | str_ends(element, "\\/"))){
        
        # Extracting first operator
        first_char <- str_sub(element, start = 1, end = 1)
        new_element <- sub('.', '', element)
        
        # Extracting last operator
        last_char <- str_sub(element, start = -1, end = -1)
        new_element <- gsub('.$', '', new_element)
        
        result <- c(result, first_char, new_element, last_char)
      }
      else if(str_starts(element, "\\*") | str_starts(element, "\\/")){
        
        # Extracting the first character of the element string
        # which should be an operator
        first_char <- str_sub(element, start = 1, end = 1)
        new_element <- sub('.', '', element)
        
        result <- c(result, first_char, new_element)
      }
      else if(str_ends(element, "\\*") | str_ends(element, "\\/")){
        
        # Extracting the last character of the element string
        # which should be an operator
        last_char <- str_sub(element, start = -1, end = -1)
        new_element <- gsub('.$', '', element)
        
        result <- c(result, new_element, last_char)
      }
      else{
        result <- c(result, element)
      }
    }
    else{
      result <- c(result, element)
    }
  }
  
  return(as.list(result))
}


# Calculating the magnitude of a given unit expression
calculate_magnitude <- function(uStr_processed){
  
  parentheses_check <- TRUE
  
  while(isTRUE(parentheses_check)){
    if(TRUE %in% str_detect(uStr_processed, "\\(")){
      parentheses_check <- TRUE
    }
    else{
      parentheses_check <- FALSE
    }
    
    if(isTRUE(parentheses_check)){
      
      position <- 0
      
      for(element in uStr_processed){
        
        position <- position + 1
        
        # Checking if element of the list is a parenthesis
        if(str_detect(element, "\\(") & str_detect(element, "\\)")){
          
          paren_pos <- position
          
          # Splitting the parenthesis based on the operators ("." or "/")
          split <- str_split(element, "(?<=[\\/\\*])|(?=[\\/\\*])", simplify = TRUE)
          
          # print("==============")
          # print(split)
          # print("==============")
          
          expression <- ""
          
          for(split_element in split){
            
            # Removing opening parenthesis
            new_element <- str_remove(split_element, "\\(")
            
            if(isTRUE(str_detect(new_element, "\\)"))){
              
              # Finding the number of closed parentheses
              n_close_parens <- str_count(new_element, "\\)")
              if(n_close_parens > 0){
                n_close_parens <- n_close_parens - 1
              }
              
              # Removing closing parenthesis
              new_element <- str_remove_all(new_element, "\\)")
            }
            
            if(!is.na(suppressWarnings(as.numeric(new_element)))){
              expression <- paste0(expression,as.numeric(new_element))
            }
            else if(new_element == "*"){
              expression <- paste0(expression,"*")
            }
            else if(new_element == "/"){
              expression <- paste0(expression,"/")
            }
            else{
          
              if(isTRUE(str_detect(new_element, "\\^"))){
                
                exponent_split <- strsplit(new_element, "\\^")
                
                position <- 0
                
                for(exp_split_element in exponent_split[[1]]){
                  position <- position + 1
                  
                  # check if exp_split_element is not a digit
                  if(is.na(suppressWarnings(as.numeric(exp_split_element)))){
                    
                    # Attempting to extract the magnitude
                    unit <- check_unit_existence(exp_split_element,
                                                 parsed_unit = TRUE)
                    exponent_split[[1]][position] <- unit$magnitude_
                  }
                }
                
                new_element <- as.character(paste(exponent_split[[1]], collapse = "^"))
                expression <- paste0(expression, new_element)
              }
              else{
                # Look for the unit inside units_df
                unit <- check_unit_existence(new_element,
                                             parsed_unit = TRUE)
                
                expression <- paste0(expression, unit$magnitude_)
              }
            }
          }
          
          # print(expression)
          
          # Parsing the expression (the result should be a digit)
          expression_result <- eval(parse(text = expression))
          
          # Adding the remaining closed parentheses to the expression result
          parentheses <- rep(")", times = n_close_parens)
          parentheses <- paste(parentheses, collapse = "")
          expression_result <- paste0(expression_result, parentheses)
          
          # print("$$$$$$$$$$$$")
          # print(expression)
          # print(expression_result)
          # print("$$$$$$$$$$$$")
          
          # Replacing the expression with the result
          uStr_processed[paren_pos] <- expression_result
        }
      }
      
      uStr_processed <- as.character(paste(uStr_processed, collapse = ""))
      uStr_processed <- processParens(uStr_processed, check_nested_parentheses(uStr_processed))
    }
    else{
      split <- str_split(uStr_processed, "(?<=[\\/\\*])|(?=[\\/\\*])", simplify = TRUE)
      
      
      # print("==============")
      # print(split)
      # print("==============")
      
      
      expression <- "" 
      
      for(split_element in split){
        
        if(isTRUE(str_detect(split_element, "\\^"))){

          exponent_split <- strsplit(split_element, "\\^")
          
          position <- 0
          
          for(exp_split_element in exponent_split[[1]]){
            position <- position + 1
            
            # check if exp_split_element is not a digit
            if(is.na(suppressWarnings(as.numeric(exp_split_element)))){
              
              # Attempting to extract the magnitude
              unit <- check_unit_existence(exp_split_element,
                                           parsed_unit = TRUE)
              exponent_split[[1]][position] <- unit$magnitude_
            }
          }
          
          new_element <- as.character(paste(exponent_split[[1]], collapse = "^"))
          expression <- paste0(expression, new_element)
        }
        else if(!is.na(suppressWarnings(as.numeric(split_element)))){
          expression <- paste0(expression, as.numeric(split_element))
        }
        else if(split_element == "*" | split_element == "/"){
          expression <- paste0(expression, split_element)
        }
        else{
          unit <- check_unit_existence(split_element,
                                       parsed_unit = TRUE)
          
          expression <- paste0(expression, unit$magnitude_)
        }
      }
      
      # print(expression)
      
      magnitude <- eval(parse(text = expression))
      return(magnitude)
    }
  }
}



# Detecting and processing parts of the string that contain
# exponents in order to indicate that they are exponents
processExponents <- function(uStr){
  
  # Extracting all strings that appear to contain an exponent
  # in either of the following formats:
  #
  # UnitDigit, Unit+Digit, Unit-Digit
  
  # pattern <- "\\b(?:[a-zA-Z]+[+|-]*\\d+)"
  pattern <- "\\b(?:[a-zA-Z]+[+|-]*\\d+)|(?:[a-zA-Z]*\\[[a-zA-Z]*\\][+|-]*\\d+)"
  exponent_matches <- unlist(str_match_all(uStr, pattern))
  
  # Extracting the remaining string parts
  remaining_string <- unlist(strsplit(uStr, pattern))
  
  # # Combining the exponent_matches with the remaining_string.
  # # The result should be exactly the same as uStr
  # result <- c(rbind(remaining_string, exponent_matches))

  # Combining the exponent_matches with the remaining_string.
  # The result should be exactly the same as uStr
  result <- character(length(uStr))
  for (i in seq_along(remaining_string)) {
    result[2*i - 1] <- remaining_string[i]
  }
  for (i in seq_along(exponent_matches)) {
    result[2*i] <- exponent_matches[i]
  }
  
  
  # Checking that the two strings are the same
  if(as.character(paste(result, collapse = "")) != uStr){
    print("String mismatch!")
    print("=== DEBUG ===")
    print(as.character(paste(result, collapse = "")))
    print(uStr)
    print(exponent_matches)
    print(remaining_string)
    print(result)
    stop("String mismatch!")
  }
  else{
    
    position <- 0
    for(element in result){
      position <- position + 1
      
      # Checking if exponent string contains a "+" (positive exponent)
      if(element %in% exponent_matches){
        if(str_detect(element, "\\+")){
          exponent_split <- strsplit(element, "\\+")
          new_element <- as.character(paste(exponent_split[[1]], collapse = "^"))
          
          result[position] <- new_element
        }
        
        # Checking if exponent string contains a "-" (negative exponent)
        else if(str_detect(element, "\\-")){
          exponent_split <- strsplit(element, "\\-")
          new_element <- as.character(paste(exponent_split[[1]], collapse = "^-"))
          
          result[position] <- new_element
        }
        
        # Assuming that exponent string doesn't contain anything (positive exponent)
        else{
          non_digits <- str_extract(element, "\\D+")
          digits <- str_extract(element, "\\d+")
          
          new_element <- paste0(non_digits, "^", digits) 
          result[position] <- new_element
        }
      }
    }
    
    return(as.character(paste(result, collapse = "")))
  }
}



# Calculating the magnitude of a given unit expression
calculate_dimension <- function(uStr_processed){
  split <- str_split(uStr_processed, "(?<=[\\/\\*])|(?=[\\/\\*])", simplify = TRUE)
  
  # print("==============")
  # print(split)
  # print("==============")
  
  expression <- "" 
  
  for(split_element in split){
    
    if(str_detect(split_element, "\\(") | str_detect(split_element, "\\)")){
      
      split_parens <- str_split(split_element, "(?<=[\\(\\)])|(?=[\\(\\)])", simplify = TRUE)
      
      parens_position <- 0
      
      for(split_parens_element in split_parens){
        
        parens_position <- parens_position + 1
        
        if(split_parens_element != "" & split_parens_element != "(" &
           split_parens_element != ")"){
          
          if(isTRUE(str_detect(split_parens_element, "\\^"))){
            
            exponent_split <- strsplit(split_parens_element, "\\^")
            
            position <- 0
            
            split_contains_unit <- FALSE
            
            for(exp_split_element in exponent_split[[1]]){
              position <- position + 1
              
              # check if exp_split_element is not a digit
              if(is.na(suppressWarnings(as.numeric(exp_split_element)))){
                
                # Attempting to extract the magnitude
                unit <- check_unit_existence(exp_split_element,
                                             parsed_unit = TRUE)
                
                exponent_split[[1]][position] <- as.character(unit$dim_[[1]])
                split_contains_unit <- TRUE
              }
              else if(!is.na(suppressWarnings(as.numeric(exp_split_element))) & !split_contains_unit){
                exponent_split[[1]][position] <- "c(0, 0, 0, 0, 0, 0, 0)"
                split_contains_unit <- FALSE
              }
            }
            
            if(split_contains_unit){
              new_element <- as.character(paste(exponent_split[[1]], collapse = "*"))
            }
            else{
              new_element <- as.character(paste(exponent_split[[1]], collapse = "+"))
            }
            
            expression <- paste0(expression, new_element)
          }
          else if(!is.na(suppressWarnings(as.numeric(split_parens_element)))){
            expression <- paste0(expression, "c(0, 0, 0, 0, 0, 0, 0)")
          }
          else if (split_parens_element == "*"){
            expression <- paste0(expression, "+")
          }
          else if (split_parens_element == "/"){
            expression <- paste0(expression, "-")
          }
          else{
            unit <- check_unit_existence(split_parens_element,
                                         parsed_unit = TRUE)
            
            expression <- paste0(expression, as.character(unit$dim_[[1]]))
          }
        }
        else{
          expression <- paste0(expression, as.character(split_parens_element))
        }
      }
    }
    else{
      if(isTRUE(str_detect(split_element, "\\^"))){
        
        exponent_split <- strsplit(split_element, "\\^")
        
        position <- 0
        
        split_contains_unit <- FALSE
        
        for(exp_split_element in exponent_split[[1]]){
          position <- position + 1
          
          # check if exp_split_element is not a digit
          if(is.na(suppressWarnings(as.numeric(exp_split_element)))){
            
            # Attempting to extract the magnitude
            unit <- check_unit_existence(exp_split_element,
                                         parsed_unit = TRUE)
            
            exponent_split[[1]][position] <- as.character(unit$dim_[[1]])
            split_contains_unit <- TRUE
          }
          else if(!is.na(suppressWarnings(as.numeric(exp_split_element))) & !split_contains_unit){
            exponent_split[[1]][position] <- "c(0, 0, 0, 0, 0, 0, 0)"
            split_contains_unit <- FALSE
          }
        }
        
        if(split_contains_unit){
          new_element <- as.character(paste(exponent_split[[1]], collapse = "*"))
        }
        else{
          new_element <- as.character(paste(exponent_split[[1]], collapse = "+"))
        }
        
        expression <- paste0(expression, new_element)
      }
      else if(!is.na(suppressWarnings(as.numeric(split_element)))){
        expression <- paste0(expression, "c(0, 0, 0, 0, 0, 0, 0)")
      }
      else if (split_element == "*"){
        expression <- paste0(expression, "+")
      }
      else if (split_element == "/"){
        expression <- paste0(expression, "-")
      }
      else{
        unit <- check_unit_existence(split_element,
                                     parsed_unit = TRUE)
        
        expression <- paste0(expression, as.character(unit$dim_[[1]]))
      }
    }
  }

  # print(expression)
  magnitude <- eval(parse(text = expression))
  # print(magnitude)
  return(magnitude)
}



# Calculating the mole exponent of a given unit expression
calculate_moleExponent <- function(uStr_processed){
  split <- str_split(uStr_processed, "(?<=[\\/\\*])|(?=[\\/\\*])", simplify = TRUE)
  
  # print("==============")
  # print(split)
  # print("==============")
  
  expression <- "" 
  
  for(split_element in split){
    
    if(str_detect(split_element, "\\(") | str_detect(split_element, "\\)")){
      
      split_parens <- str_split(split_element, "(?<=[\\(\\)])|(?=[\\(\\)])", simplify = TRUE)
      
      parens_position <- 0
      
      for(split_parens_element in split_parens){
        
        parens_position <- parens_position + 1
        
        if(split_parens_element != "" & split_parens_element != "(" &
           split_parens_element != ")"){
          
          if(isTRUE(str_detect(split_parens_element, "\\^"))){
            
            exponent_split <- strsplit(split_parens_element, "\\^")
            
            position <- 0
            
            split_contains_unit <- FALSE
            
            for(exp_split_element in exponent_split[[1]]){
              position <- position + 1
              
              # check if exp_split_element is not a digit
              if(is.na(suppressWarnings(as.numeric(exp_split_element)))){
                
                # Attempting to extract the mole exponent
                unit <- check_unit_existence(exp_split_element,
                                             parsed_unit = TRUE)
                
                exponent_split[[1]][position] <- as.character(unit$moleExp_[[1]])
                split_contains_unit <- TRUE
              }
              else if(!is.na(suppressWarnings(as.numeric(exp_split_element))) & !split_contains_unit){
                exponent_split[[1]][position] <- "0"
                split_contains_unit <- FALSE
              }
              else{
                exponent_split[[1]][position] <- "1"
                split_contains_unit <- TRUE
              }
            }
            
            if(split_contains_unit){
              new_element <- as.character(paste(exponent_split[[1]], collapse = "*"))
            }
            else{
              new_element <- as.character(paste(exponent_split[[1]], collapse = "+"))
            }
            
            expression <- paste0(expression, new_element)
          }
          else if(!is.na(suppressWarnings(as.numeric(split_parens_element)))){
            expression <- paste0(expression, "0")
          }
          else if (split_parens_element == "*"){
            expression <- paste0(expression, "+")
          }
          else if (split_parens_element == "/"){
            expression <- paste0(expression, "-")
          }
          else{
            unit <- check_unit_existence(split_parens_element,
                                         parsed_unit = TRUE)
            
            expression <- paste0(expression, as.character(unit$moleExp_[[1]]))
          }
        }
        else{
          expression <- paste0(expression, as.character(split_parens_element))
        }
      }
    }
    else{
      if(isTRUE(str_detect(split_element, "\\^"))){
        
        exponent_split <- strsplit(split_element, "\\^")
        
        # print(exponent_split)
        
        position <- 0
        
        split_contains_unit <- FALSE
        
        for(exp_split_element in exponent_split[[1]]){
          position <- position + 1
          
          # check if exp_split_element is not a digit
          if(is.na(suppressWarnings(as.numeric(exp_split_element)))){
            
            # Attempting to extract the mole exponent
            unit <- check_unit_existence(exp_split_element,
                                         parsed_unit = TRUE)
            
            exponent_split[[1]][position] <- as.character(unit$moleExp_[[1]])
            split_contains_unit <- TRUE
          }
          else if(!is.na(suppressWarnings(as.numeric(exp_split_element))) & !split_contains_unit){
            exponent_split[[1]][position] <- "0"
            split_contains_unit <- FALSE
          }
          else{
            exponent_split[[1]][position] <- "1"
            split_contains_unit <- TRUE
          }
        }
        
        if(split_contains_unit){
          new_element <- as.character(paste(exponent_split[[1]], collapse = "*"))
        }
        else{
          new_element <- as.character(paste(exponent_split[[1]], collapse = "+"))
        }
        
        expression <- paste0(expression, new_element)
      }
      else if(!is.na(suppressWarnings(as.numeric(split_element)))){
        expression <- paste0(expression, "0")
      }
      else if (split_element == "*"){
        expression <- paste0(expression, "+")
      }
      else if (split_element == "/"){
        expression <- paste0(expression, "-")
      }
      else{
        unit <- check_unit_existence(split_element,
                                     parsed_unit = TRUE)
        
        expression <- paste0(expression, as.character(unit$moleExp_[[1]]))
      }
    }
  }
  
  # print(expression)
  moleExp <- eval(parse(text = expression))
  # print(moleExp)
  return(moleExp)
}



# Identifying whether the unit expression contains non-ratio units (cnv_)
# and also checking for the expression's validity

identify_cnv <- function(uStr_processed){
  
  # Setting the unit counters
  n_units <- 0
  n_cnv_units <- 0
  
  # Setting the default cnv_ to NA
  cnv_ <- NA
  
  csCode_ <- NA
  
  # Detecting operators inside the unit expression
  operators <- str_detect(uStr_processed, "\\/|\\*")
  division <- str_detect(uStr_processed, "\\/")
  
  # Splitting expression based on the operators (/,*)
  split <- str_split(uStr_processed, "(?<=[\\/\\*])|(?=[\\/\\*])", simplify = TRUE)
  
  # print("==============")
  # print(split)
  # print("==============")
  
  
  for(split_element in split){
    
    # Checking if the split element contains an open or closed parenthesis
    if(str_detect(split_element, "\\(") | str_detect(split_element, "\\)")){
      
      # Splitting on the parenthesis level
      split_parens <- str_split(split_element, "(?<=[\\(\\)])|(?=[\\(\\)])", simplify = TRUE)
      
      # print("==============")
      # print(split_parens)
      # print("==============")
      
      for(split_parens_element in split_parens){
        
        # Checking the split elements that do not contain a whitespace
        # or a parenthesis
        if(split_parens_element != "" & split_parens_element != "(" &
           split_parens_element != ")"){
          
          # Checking if parenthesis element contains an exponent
          if(isTRUE(str_detect(split_parens_element, "\\^"))){
            
            # Splitting based on the exponent
            exponent_split <- strsplit(split_parens_element, "\\^")
            
            # print(exponent_split)
            
            for(exp_split_element in exponent_split[[1]]){
              
              # Checking if the exponent element is non-numeric
              if(is.na(suppressWarnings(as.numeric(exp_split_element)))){
                
                # Checking if the exponent element is a unit that exists
                # inside the cnv_units.
                # If it exists, increase the n_cnv_units counter
                if(exp_split_element %in% cnv_units$csCode_ |
                   exp_split_element %in% cnv_units$ciCode_){
                  n_cnv_units <- n_cnv_units + 1
                  
                  # Filtering the cnv_units with the found unit to extract
                  # its cnv_ value
                  cnv_unit <- cnv_units %>% filter(csCode_ == exp_split_element)
                  if(nrow(cnv_unit) < 1){
                    cnv_unit <- cnv_units %>% filter(ciCode_ == exp_split_element)
                  }
                  cnv_ <- cnv_unit$cnv_
                  csCode_ <- cnv_unit$csCode_
                }
                # If it doesn't exist, check whether it exists inside the
                # units_df as a unit
                # If it exists, increase the n_units counter
                else{
                  unit <- check_unit_existence(exp_split_element,
                                               parsed_unit = TRUE)
                  if(nrow(unit) == 1){
                    n_units <- n_units + 1
                  }
                  # if it doesn't exist, throw an error message
                  else{
                    stop(paste0("ERROR UNIT NOT FOUND: ",exp_split_element))
                  }
                }
              }
            }
          }
          # If parenthesis element doesn't include an exponent,
          # check if the parenthesis element is a unit that exists
          # inside the cnv_units
          else{
            # If it exists, increase the n_cnv_units counter
            if(is.na(suppressWarnings(as.numeric(split_parens_element)))){
              if(split_parens_element %in% cnv_units$csCode_ |
                 split_parens_element %in% cnv_units$ciCode_){
                n_cnv_units <- n_cnv_units + 1
                
                # Filtering the cnv_units with the found unit to extract
                # its cnv_ value
                cnv_unit <- cnv_units %>% filter(csCode_ == split_parens_element)
                if(nrow(cnv_unit) < 1){
                  cnv_unit <- cnv_units %>% filter(ciCode_ == split_parens_element)
                }
                cnv_ <- cnv_unit$cnv_
                csCode_ <- cnv_unit$csCode_
              }
              # If it doesn't exist, check whether it exists inside the
              # units_df as a unit
              # If it exists, increase the n_units counter
              else{
                unit <- check_unit_existence(split_parens_element,
                                             parsed_unit = TRUE)
                if(nrow(unit) == 1){
                  n_units <- n_units + 1
                }
                # if it doesn't exist, throw an error message
                else{
                  stop(paste0("ERROR UNIT NOT FOUND: ",split_parens_element))
                }
              }
            }
          }
        }
      }
    }
    # If split element doesn't contain a parenthesis
    else{
      # Check if split element contains an exponent
      if(isTRUE(str_detect(split_element, "\\^"))){
        
        # Splitting based on the exponent
        exponent_split <- strsplit(split_element, "\\^")
        
        # print(exponent_split)
        
        for(exp_split_element in exponent_split[[1]]){
          
          # Checking if the exponent element is non-numeric
          if(is.na(suppressWarnings(as.numeric(exp_split_element)))){
            
            # Checking if the exponent element is a unit that exists
            # inside the cnv_units.
            # If it exists, increase the n_cnv_units counter
            if(exp_split_element %in% cnv_units$csCode_ |
               exp_split_element %in% cnv_units$ciCode_){
              n_cnv_units <- n_cnv_units + 1
              
              # Filtering the cnv_units with the found unit to extract
              # its cnv_ value
              cnv_unit <- cnv_units %>% filter(csCode_ == exp_split_element)
              if(nrow(cnv_unit) < 1){
                cnv_unit <- cnv_units %>% filter(ciCode_ == exp_split_element)
              }
              cnv_ <- cnv_unit$cnv_
              csCode_ <- cnv_unit$csCode_
            }
            # If it doesn't exist, check whether it exists inside the
            # units_df as a unit
            # If it exists, increase the n_units counter
            else{
              unit <- check_unit_existence(exp_split_element,
                                           parsed_unit = TRUE)
              if(nrow(unit) == 1){
                n_units <- n_units + 1
              }
              # if it doesn't exist, throw an error message
              else{
                stop(paste0("ERROR UNIT NOT FOUND: ",exp_split_element))
              }
            }
          }
        }
      }
      # If split element doesn't include an exponent,
      # check if the split element is a unit that exists
      # inside the cnv_units
      else{
        if(is.na(suppressWarnings(as.numeric(split_element)))){
          
          # If it exists, increase the n_cnv_units counter
          if(split_element %in% cnv_units$csCode_ |
             split_element %in% cnv_units$ciCode_){
            n_cnv_units <- n_cnv_units + 1
            
            # Filtering the cnv_units with the found unit to extract
            # its cnv_ value
            cnv_unit <- cnv_units %>% filter(csCode_ == split_element)
            if(nrow(cnv_unit) < 1){
              cnv_unit <- cnv_units %>% filter(ciCode_ == split_element)
            }
            cnv_ <- cnv_unit$cnv_
            csCode_ <- cnv_unit$csCode_
          }
          
          # If it doesn't exist, check that split element is not an operator,
          # then check whether it exists inside the units_df as a unit
          else{
            if(split_element != "/" & split_element != "*"){
              unit <- check_unit_existence(split_element,
                                           parsed_unit = TRUE)
              
              # If it exists, increase the n_units counter
              if(nrow(unit) == 1){
                n_units <- n_units + 1
              }
              # if it doesn't exist, throw an error message
              else{
                stop(paste0("ERROR UNIT NOT FOUND: ",split_element))
              }
            }
          }
        }
      }
    }
  }
  
  # Throw an error if there are multiple cnv_ units inside the expression
  if(n_cnv_units > 1){
    stop("Invalid Unit Expression")
  }
  
  # Throw an error if there are both cnv_ and other units inside the expression
  if(n_cnv_units == 1 & n_units > 0 & isTRUE(operators)){
    stop("Unable to multiply or divide with non-ratio units")
  }
  
  # Throw an error if a division is occurring inside the expression 
  if(n_cnv_units == 1 & isTRUE(division)){
    stop("Attempt to divide with a non-ratio unit")
  }
  
  # Return TRUE only if the expression contains a single cnv_ unit
  if(n_cnv_units == 1 & n_units == 0){
    return(list(cnv_exists = TRUE,
                cs_code = csCode_,
                cnv_value = cnv_))
  }
  else{
    return(list(cnv_exists = FALSE,
                cs_code = csCode_,
                cnv_value = cnv_))
  }
}





# Calculating the magnitude_ and cnvPfx_ of cnv unit expressions
calculate_cnv_mag_pfx <- function(uStr, cnv_csCode_){
  
  # Formatting cnv_csCode_ by revealing the escape characters
  # "\\" in order to properly process cases where 
  # cnv_csCode_ contains "[...]"
  format_cnv_csCode_ <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", cnv_csCode_, perl = TRUE)
  
  # Establishing the pattern used for extracting the cnv_ from
  # the uStr expression containg any potential exponents
  pattern <- paste0(format_cnv_csCode_,"\\^?[\\+|\\-]?\\d*")
  
  # Extracting the cnv_ with the exponents from uStr
  extract <- str_extract(uStr, pattern)
  
  # Filtering the cnv_units for the cnv used inside the expression
  # and extracting its magnitude and cnvPfx_ value (a unit should ALWAYS be found)
  unit <- cnv_units %>% filter(csCode_ == cnv_csCode_)
  cnv_unit_magnitude <- unit$magnitude_
  cnvPfx_ <- unit$cnvPfx_
  
  # Replacing the cnv_ with it's magnitude and parsing the expression
  new_extract <- str_replace(extract, format_cnv_csCode_, as.character(cnv_unit_magnitude))
  cnv_magnitude <- eval(parse(text = new_extract))
  
  
  
  # Splitting the uStr expression based on the defined pattern
  substring <- str_split(uStr, pattern)
  
  expression <- ""
  
  for(string in substring[[1]]){
    # Removing any opening parenthesis at the end of the string
    if(str_ends(string, "\\(")){
      new_string <- str_replace(string, "\\(", "")
      expression <- paste0(expression, new_string)
    }
    # Removing any closing parentheses at the beginning of the string,
    # as well as any operator
    else if(str_starts(string, "\\)?[\\*|\\/]*")){
      new_string <- str_replace(string, "\\)?[\\*|\\/]*", "")
      expression <- paste0(expression, new_string)
    }
    else{
      expression <- paste0(expression, string)
    }
  }
  
  expression <- paste0(expression,"*",cnvPfx_)
  
  # print(uStr)
  # print(extract)
  # print(pattern)
  # print(substring)
  # print(expression)
  
  cnvPfx_ <- eval(parse(text = expression)) 
  
  return(list(cnv_magnitude = cnv_magnitude,
              cnvPfx_ = cnvPfx_))
}







