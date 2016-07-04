# GTEx eQTL database loading script
#
# PLAN:  
#    1. for each eQTL file (each a single tissue), extract the tissue name
#    2. load eQTL data into staging table  
#    3. update the eQTL table. This requires a (large) join to both dimTissue and dimGene.
#
# NOTES:
# ------
#   - ASSUMES that the data warehouse is already populated with dimTissue, dimGene, dimDataSource and factExpression.
#   - user should enter appropriate username, password and hostname 
#     (system calls to db_query() requires blank password at this stage
#      WARNING: remember to add password to user afterwards)
#   - this ETL script should be saved in the parent directory above where the 
#     raw eQTL files are saved. (i.e. on the database server!)  
#   - may be run either locally, or on the database server's host  
#   - originally queries were passed to the db via RMySQL::dbGetQuery(), however this 
#     caused excessive row locks. Seems to work when making system calls (see db_query() below)
#
# Nick Burns
# June 2016

library(RMySQL)
setwd("./eQTL")

# GLOBAL VARIABLES
# SQL Connection settings - you need to put appropriate values in here.
mysql_user <- ""    
mysql_password <- ""
mysql_host <- ""

eqtl_file_pattern <- "_Analysis_cis-eQTLs.txt"
eqtl_files <- list.files(".", pattern = eqtl_file_pattern)

# HELPER FUNCTIONS
strip_suffix <- function (x, suffix = "_Analysis_cis-eQTLs.txt") {
    
    # Strips the end off a file name (or string)
    tissue <- substr(x, start = 1, stop = nchar(x) - nchar(suffix))
    return (tissue)
}

reset_staging <- function (conn) {
    query <- "call qtl_reset_staging();"
}
etl <- function (eqtl_file) {
    
    # Performs a bulk load, ignoring the header in each file
    query <- sprintf("LOAD DATA LOCAL INFILE '%s'
                     INTO TABLE eqtl_staging
                     IGNORE 1 LINES;", eqtl_file)
    return (query)
}
pop_fact <- function (tissue_id, source_id) {
    query <- sprintf("
                     LOCK TABLES factQTL WRITE, eqtl_staging WRITE, eqtl_staging as stage WRITE, dimGene WRITE, dimGene as gene WRITE;
                     CALL qtl_populate_fact(%s, %s);
                     UNLOCK TABLES;", tissue_id, source_id)
}
db_query <- function (query, username = "nickburns", host = "biocvisg0.otago.ac.nz", db = "eQTL_dw") {
    execute <- sprintf('mysql -u %s -h %s -D %s -e "%s"',
                       username, host, db, query)
    return(execute)
}
db_commit <- function () {
    system(db_query("commit;"))
}
toggle_config <- function (config_setting, state = 1) {
    # toggles the state of a given configuration setting
    # NOTE: default state is to ENABLE a setting
    query <- sprintf("SET %s = %s;", config_setting, state)
    system(db_query(query))
}

# main loop...
master_start <- Sys.time()

# BULK INSERT optimisation setting
toggle_config("autocommit", state = 0)
toggle_config("foreign_key_checks", state = 0)
toggle_config("sql_log_bin", state = 0)
toggle_config("unique_checks", state = 0)

for (eQTL_file in eqtl_files) {
    
    print(sprintf("-------------    %s    -------------", eQTL_file))
    print("")
    lcl_start <- Sys.time()
    
    # extract tissue name
    tissue <- strip_suffix(eQTL_file)

    # load eQTL data in staging table
    print("    ... bulk load into staging")
    system(db_query(reset_staging()))
    system(db_query(etl(eQTL_file)))
    db_commit()
    
    # populate factQTL...
    
    
    # informational only
    lcl_end <- Sys.time()
    print(sprintf("Time for %s: %s", eQTL_file, lcl_end - lcl_start))
    print("")
    print("")
}
# finally, truncate the staging table
system(db_query(reset_staging()))

# reset BULK INSERT optimisation settings
toggle_config("autocommit")
toggle_config("foreign_key_checks")
toggle_config("sql_log_bin")
toggle_config("unique_checks")

master_end <- Sys.time()
print(sprintf("Total time:  %s", master_end - master_start))
