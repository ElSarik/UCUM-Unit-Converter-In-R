# UCUM Unit Converter In R
The repository contains an implementation of the [UCUM-LHC unit converter](https://ucum.nlm.nih.gov/ucum-lhc/demo.html) in R.

The converter was designed as part of an upcoming _**(as of 12 March 2024)**_ R package called "lab2clean", which provides automation and standardization of cleaning restrospective clinical lab data. The provided code is a standalone version of the conversion functions used in the R package.


## Disclaimer
The implementation is NOT meant to be a direct translation of the original [UCUM-LHC source code](https://github.com/lhncbc/ucum-lhc) from Javascript to R.

A reverse engineering approach was used in order to understand the tool's functionality, which led to simplification of the original code, while still maintaining the core features of the tool.

The conversion results were validated with the original code and tool, but due to the complexity of the UCUM unit expressions there may be some conversion mistakes which have not yet been identified.


## Features
* Conversion between any **valid** UCUM units, including units which require a moleqular weight
* Custom made unit parser for performing conversions between complex unit expressions 


## Requirements
The code was developed and tested using the following:
* R version 4.3.0 (2023-04-21 ucrt) -- "Already Tomorrow"
* tidyverse Version 2.0.0


## Usage
Run the converter from within run.R
