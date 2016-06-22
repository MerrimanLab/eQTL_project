# GTEx eQTL database loading script
#
# PLAN:  
#    1. for each eQTL file (each a single tissue), extract the tissue name
#    2. load eQTL data into staging table  
#    3. update tissue dimension table
#    4. update gene dimension table  (already loaded)
#    5. update the eQTL table
#
# NOTES:
# ------
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
mysql_user <- "nickburns"    
mysql_password <- ""
mysql_host <- "biocvisg0.otago.ac.nz"

eqtl_file_pattern <- "_Analysis_cis-eQTLs.txt"
eqtl_files <- list.files(".", pattern = eqtl_file_pattern)

# HELPER FUNCTIONS
strip_suffix <- function (x, suffix = "_Analysis_cis-eQTLs.txt") {
    
    # Strips the end off a file name (or string)
    tissue <- substr(x, start = 1, stop = nchar(x) - nchar(suffix))
    return (tissue)
}
pop_tissue_dimension <- function (x) {
    query <- sprintf("INSERT INTO dimTissue VALUES (DEFAULT, '%s')", x)
    return (query)
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
pop_gene_dimension <- function () {
    query <- "
                LOCK TABLES dimGene WRITE, eQTL_staging WRITE, dimGene as g WRITE, eQTL_staging as stage WRITE ;
                INSERT INTO dimGene (ensembl_id, chromosome)
                SELECT DISTINCT stage.ensembl_id, stage.chromosome
                FROM (
                        SELECT 
                            substring_index(ensembl_id, '.', 1) as ensembl_id,
                            substring_index(snp_id,'_',1) as chromosome
                        FROM eQTL_staging
                     ) as stage
                LEFT OUTER JOIN dimGene g on stage.ensembl_id = g.ensembl_id
                WHERE g.ensembl_id IS NULL;
                UNLOCK TABLES;
    "
    return (query)
}
get_tissue_id <- function (tissue) {
    query <- sprintf(
        "SELECT * from dimTissue
        WHERE tissue_description = '%s';", tissue)
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
    query <- sprintf("SET %s %s;", config_setting, state)
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

    # update tissue dim
    print("    ... updating dimTissue")
    system(db_query(pop_tissue_dimension(tissue)))
    db_commit()
    
    # load eQTL data in staging table
    print("    ... bulk load into staging")
    system(db_query(reset_staging()))
    db_commit()

    system(db_query(etl(eQTL_file)))
    db_commit()
    
    system(db_query("create index idx_staging on eqtl_staging (pvalue) using hash"))
    db_commit()
    
#     # update gene dim
#     print("    ... updating dimGene")
#     system(db_query(pop_gene_dimension()))
#     # NOTE: this has already been bulk loaded.
    
    # update factQTL
    print("    ... updating fact table")
    conn <- RMySQL::dbConnect(RMySQL::MySQL(),
                              username = mysql_user,
                              password = mysql_password,
                              host = mysql_host,
                              dbname = "eQTL_dw")
    tissue_id <- dbGetQuery(conn, get_tissue_id(tissue))
    dbDisconnect(conn)
    
    system(db_query(pop_fact(tissue_id$tissue_id, 1)))
    db_commit()
    
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
