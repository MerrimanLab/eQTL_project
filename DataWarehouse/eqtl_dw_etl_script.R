# GTEx eQTL database loading script
#
# PLAN:  
#    1. for each eQTL file (each a single tissue), extract the tissue name
#    2. load eQTL data into staging table  
#    3. update tissue dimension table
#    4. update gene dimension table  
#    5. update the eQTL table
#
# NOTES:
# ------
#   - user should enter appropriate username, password and hostname 
#     (system calls to db_query() requires blank password at this stage
#      WARNING: remember to add password to user afterwards)
#   - this ETL script should be saved in the parent directory above where the 
#     raw eQTL files are saved.  
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
pop_tissue_dimension <- function (x) {
    query <- sprintf("INSERT INTO dimTissue VALUES (DEFAULT, '%s')", x)
    return (query)
}
reset_staging <- function (conn) {
    dbGetQuery(conn, "delete from eQTL_staging;")
}
etl <- function (eqtl_file) {
    
    # Performs a bulk load, ignoring the header in each file
    query <- sprintf("LOAD DATA LOCAL INFILE '%s'
                     INTO TABLE eQTL_staging
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
                WHERE g.ensembl_id IS NULL limit 500;
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
pop_fact <- function (tissue_id) {
    query <- sprintf("
    LOCK TABLES eQTL_staging WRITE, factQTL WRITE;
    INSERT INTO factQTL (ensembl_id, tissue, chromosome, build_37_pos, beta, tstat, pvalue, source_name)
    SELECT 
        substring_index(ensembl_id, '.', 1),
        '%s',
        substring_index(snp_id, '_', 1),
        substring_index(substring_index(snp_id, '_', 2), '_', -1),
        beta,
        tstat,
        pvalue,
        'GTEx'
    FROM eQTL_staging;
    UNLOCK TABLES;
    ", tissue_id)
}
db_query <- function (query, username = "nickburns", host = "biocvisg0.otago.ac.nz", db = "eQTL_dw") {
    execute <- sprintf('mysql -u %s -h %s -D %s -e "%s"',
                       username, host, db, query)
    return(execute)
}

# main loop...
master_start <- Sys.time()
for (eQTL_file in eqtl_files) {
    
    print(sprintf("-------------    %s    -------------", eQTL_file))
    print("")
    lcl_start <- Sys.time()
    
    # extract tissue name
    tissue <- strip_suffix(eQTL_file)

    # update tissue dim
    print("    ... updating dimTissue")
    system(db_query(pop_tissue_dimension(tissue)))
    
    # load eQTL data in staging table
    print("    ... bulk load into staging")
    system(db_query("delete from eQTL_staging;"))
    system(db_query(etl(eQTL_file)))
    
    # update gene dim
    print("    ... updating dimGene")
    system(db_query(pop_gene_dimension()))
    
    # update factQTL
    print("    ... updating fact table")
    conn <- RMySQL::dbConnect(RMySQL::MySQL(),
                              username = mysql_user,
                              password = mysql_password,
                              host = mysql_host,
                              dbname = "eQTL_dw")
    tissue_id <- dbGetQuery(conn, get_tissue_id(tissue))
    dbDisconnect(conn)
    system(db_query(pop_fact(tissue_id$tissue_id)))
    
    lcl_end <- Sys.time()
    print(sprintf("Time for %s: %s", eQTL_file, lcl_end - lcl_start))
    print("")
    print("")
    
}

master_end <- Sys.time()
print(sprintf("Total time:  %s", master_end - master_start))
