# GTEx eQTL database loading script
#
# PLAN:  
#    1. for each eQTL file (each a single tissue), extract the tissue name
#    2. load eQTL data into staging table  
#    3. update tissue dimension table
#    4. update gene dimension table  
#    5. update the eQTL table
#
# The process is scripted here, largely just so I can loop and extract tissue
# names easily. However, all heavy lifting is handled by calls to mysql.
#
# Nick Burns
# June 2016

library(RMySQL)
setwd("/mnt/DataDrive/gEXPR_eQTL_Datasets/GTEXData/eQTL/")

# GLOBAL VARIABLES
# SQL Connection settings - you need to put appropriate values in here.
mysql_user <- ""    
mysql_password <- ""
mysql_host <- ""

eqtl_file_pattern <- "_Analysis_cis-eQTLs.txt.gz"
eqtl_files <- list.files(".", pattern = eqtl_file_pattern)

# HELPER FUNCTIONS
strip_suffix <- function (x, suffix = "_Analysis_cis-eQTLs.txt.gz") {
    
    # Strips the end off a file name (or string)
    tissue <- substr(x, start = 1, stop = nchar(x) - nchar(suffix))
    return (tissue)
}
pop_tissue_dimension <- function (x) {
    query <- sprintf("INSERT INTO dimTissue VALUES (DEFAULT, 'x')")
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
    query <- "INSERT INTO dimGene (ensembl_id, chromosome)
                SELECT DISTINCT stage.ensembl_id, stage.chromosome
                FROM (
                        SELECT 
                            substring_index(ensembl_id, '.', 1) as ensembl_id,
                            substring_index(snp_id,'_',1) as chromosome
                        FROM eQTL_staging
                     ) as stage
                LEFT OUTER JOIN dimGene g on stage.ensembl_id = g.ensembl_id
                WHERE g.ensembl_id IS NULL;
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
    ", tissue_id)
}


# main loop...
master_start <- Sys.time()
conn <- RMySQL::dbConnect(RMySQL::MySQL(),
                          username = mysql_user,
                          password = mysql_password,
                          host = mysql_host,
                          dbname = "eQTL_dw")
for (eQTL_file in eqtl_files) {
    
    print(sprintf("-------------    %s    -------------", eQTL_file))
    print("")
    lcl_start <- Sys.time()
    
    # lcl_eqtl_file
    lcl_file <- strip_suffix(eQTL_file, suffix = ".gz")
    
    # extract tissue name
    tissue <- strip_suffix(eQTL_file)
    
    # unzip the file
    if (!file.exists(lcl_file)) {
        print("... unzipping ...")
        system(sprintf("gunzip %s", eQTL_file))
    }
    # update tissue dim
    dbGetQuery(conn, pop_tissue_dimension(tissue))
    
    # load eQTL data in staging table
    reset_staging(conn)
    dbGetQuery(conn, etl(lcl_file))
    
    # update gene dim
    dbGetQuery(conn, pop_gene_dimension())
    
    # update factQTL
    tissue_id <- dbGetQuery(conn, get_tissue_id(tissue))
    dbGetQuery(conn, pop_fact(tissue_id))
    
    lcl_end <- Sys.time()
    print(sprintf("Time for %s: %s", eQTL_file, lcl_end - lcl_start))
    print("")
    print("")
    
}
dbDisconnect(conn)

master_end <- sys.time()
print(sprintf("Total time:  %s"), master_end - master_start)
