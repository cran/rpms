###############################################################################
#' Recursive Partitioning for Modeling Survey Data (rpms)
#'
#' This package provides a function \code{rpms} to produce an \code{rpms} object 
#' and method functions that operate on them. 
#' The \code{rpms} object is a representation of a regression tree achieved
#' by recursively partioning the dataset, fitting the specified linear model
#' on each node seperately.
#' The recursive partitioning algorithm has an unbiased variable selection
#' and accounts for the sample design.
#' The algorithm accounts for one-stage of stratification and clustering as
#' well as unequal probability of selection.
#' This version does not handle missing values, so only complete cases of a 
#' dataset is used.
#' 
#' @docType package
#' @name rpms
#' 
#' @useDynLib rpms
#' @importFrom Rcpp sourceCpp
NULL

#############################################################################
#' CE Consumer expenditure data (first quarter of 2014)
#'
#' A dataset containing consumer unit characteristics, assets and expenditure 
#' data from the Bureau of Labor Statistics' Consumer Expenditure Survey 
#' public use interview data file.
#' 
#' @format A data frame with 6483 rows and 61 variables:
#' @section Location and sample-design variables:
#' \describe{
#' \item{NEWID}{Consumer unit identifying variable}
#' \item{PSU}{Primary Sampling Unit code}
#'  \item{CID}{Cluster Identifier for all clusters, (defined using PSU, 
#'            REGION, STATE, and POPSIZE)}
#' \item{FINLWT21}{Final sample weight to make inference to total population}
#' \item{POPSIZE}{Population size of PSU 1-biggest 5-smallest}
#' \item{REGION}{Region code: 1 Northeast; 2 Midwest; 3 South; 4 West}
#' \item{STATE}{State FIPS code}
#' \item{BLS_URBN}{Urban = 1, Rural = 2}
#' \item{SMSASTAT}{Is CU in a MSA: 1 Yes; 2 No}}
#' 
#' @section Household variables:
#' \describe{
#' \item{HH_CU_Q}{Number of CU in household}
#' \item{FAM_TYPE}{CU code based on relationship of members to reference person  
#'          (children incldue blood-related, step and adopted):
#'          1 Married Couple only; 
#'          2 Married Couple, children (oldest < 6 years old);
#'          3 Married Couple, children (oldest 6 to 17 years old);
#'          4 Married Couple, children (oldest > 17 years old);
#'          5 All other Married Couple CUs
#'          6 One parent (male), children (at least one child < 18 years old);
#'          7 One parent (female), children (at least one child < 18 years old);
#'          8 Single consumers; 9 Other CUs}
#' \item{FAM_SIZE}{Number of members in CU}
#' \item{AS_COMP1}{Number of males >16 yrs old}
#' \item{AS_COMP2}{Number of females >16 yrs old}
#' \item{PERSLT18}{Number of people <18 yrs old}
#' \item{PERSOT64}{Number of people >64 yrs old}
#' \item{NO_EARNR}{Number of earners}}
#' 
#' @section Housing and transportation:
#' \describe{
#' \item{CUTENURE}{Housing tenure: 1 Owned with mortgage; 
#'          2 Owned without mortgage
#'          3 Owned mortgage not reported; 4 Rented; 
#'          5 Occupied without payment of cash rent; 6 Student housing}
#' \item{ROOMSQ}{Number of rooms, including finished living areas 
#'          and excluding all baths}
#' \item{BATHRMQ}{Number of complete bathrooms in the unit}
#' \item{BEDROOMQ}{Tota number of bedrooms in the unit}         
#' \item{VEHQ}{Number of owned vehicles}
#' \item{VEHQL}{Number of leased vehicles}}
#' 
#' @section Reference person:
#' \describe{
#' \item{AGE_REF}{Age of reference person}
#' \item{EDUC_REF}{Education level of reference person coded: 
#'                  00 None; 10 1st-8th Grade; 11 some HS; 12 HS; 
#'                  13 Some college; 14 AA degree; 15 Bachelors degree; 
#'                  16 Advanced degree}
#' \item{REF_RACE}{Race code of reference person: 1 White; 2 Black; 
#'              3 Native American; 4 Asian; 5 Pacific Islander; 6 Multi-race}
#' \item{SEX_REF}{Male = 1; Female = 2}
#' \item{HISP_REF}{Hispanic = 1; Not Hispanic = 2}
#' \item{HIGH_EDU}{Highest level education level attained by anyone within CU: 
#'                  00 None; 10 1st-8th Grade; 11 some HS; 12 HS; 
#'                  13 Some college; 14 AA degree; 15 Bachelors degree; 
#'                  16 Advanced degree}}
#' 
#' @section Labor status variables:             
#' \describe{     
#' \item{INC_HRS1}{Number of hours usually worked per week by reference person}
#' \item{INCNONW1}{Reason reference person did not work during the past 12 
#'                  months:
#'                  1 Retired; 2 Home maker; 3 School; 4 health;
#'                  5 Unable to find work; 6 Doing something else}}
#' 
#' @section Income variables:
#' \describe{                  
#' \item{FINCBTAX}{Amount of CU income before taxes in past 12 months}
#' \item{FINCBTAX_I}{Imputation indicator for FINCBTAX: 0-reported; 1-imputed}
#' \item{FINCBT_X}{Flag for FINCBTAX: D-reported value; T-top coded} 
#' \item{FINCATAX}{Amount of CU income after taxes in past 12 months}
#' \item{FINCAT_X}{Flag for FINCATAX: D-reported value; T-top coded}   
#' \item{FSALARYX}{Amount of wage and salary income, before deductions, 
#'                received by all CU members in past 12 months}
#' \item{FSALARY_I}{Imputation indicator for FSALARYX: 0-reported; 1-imputed}           
#' \item{FSAL_RYX}{Flag for FSAL_RYX: D-reported; T-top coded}
#' \item{FRRETIRX}{Amount of Social Security and Railroad Retirement income}
#' \item{FRRET_I}{Imputation indicator for FRRETIRX: 0-reported; 1-imputed}
#' \item{WELFAREM}{Total amount of income from public assistance or welfare, 
#'                   including money from job training grants, received by ALL CU 
#'                   members during the past 12 months}
#' \item{WELFARE_I}{Imputation indicator for WELFAREM: 0-reported; 1-imputed}
#' \item{NETRENTX}{Amount of net rental income}
#' \item{NETRENT_I}{Imputation indicator for NETRENTX: 0-reported; 1-imputed}
#' \item{NETR_NTX}{Flag for NETRENTX: A-valid blank; C-refusal; 
#'                 D-reported value; T-value is top coded}
#' \item{ROYESTX}{Amount income received in royalties, estates and trusts}
#' \item{ROYESTX_I}{Imputation indicator for ROYESTX: 0-reported; 1-imputed}               
#' \item{ROYESTX_}{Flag for ROYESTX: A-valid blank; C-refusal; D-reported value;
#'                 T-value is top coded}}
#' @section Assets: 
#' \describe{
#' \item{IRAYRX}{value of all retirement accounts one year ago today}
#' \item{IRAYRX_}{Flag for IRAYRX: A-valid blank; C-refusal; D-reported value;
#'                 T-value is top coded}
#' \item{LIQUDYRX}{Value of all checking, savings, money market accounts, and
#'                 certificates of deposit one year ago today}
#' \item{LIQU_YRX}{Flag for LIQUDYRX: A-valid blank; C-refusal; D-reported value;
#'                 T-value is top coded}
#' \item{STOCKB}{Range of total value of all directly-held stocks, bonds, 
#'              and mutual funds: 1 $0 - $1999; 2 $2,000 - $9,999; 
#'              3 $10,000 - $49,999; 4 $50,000 - $199,999; 
#'              5 $200,000 - $449,999; 6 $450,000 and over}}
#' @section Expenditures:
#' \describe{            
#' \item{FINDRETX}{Amount of money put in an individual retirement plan, such 
#'                  as an IRA or Keogh, by all CU members in past 12 months}
#' \item{FIND_ETX}{Flag for FINDRETX: D-reported; T-top coded}                 
#' \item{TOTXEST}{Estimated total taxes paid}
#' \item{FDHOMECQ}{Expenditure on food at home this quarter}
#' \item{FDAWAYCQ}{Expenditure on food away from home this quarter}
#' \item{ALCBEVCQ}{Expenditure on alcoholic beverages this quarter}
#' \item{GASMOCQ}{Expenditure on gasoline and motor oil this quarter}
#' } #end describe
#' @source \url{http://www.bls.gov/cex/pumd_data.htm}
#' 
#' @seealso For more information see 
#' \url{http://www.bls.gov/cex/2015/csxintvw.pdf}
"CE"