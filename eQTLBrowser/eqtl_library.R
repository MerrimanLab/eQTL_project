# eQTL Browser
# eqtl_library.R
#
# Logic code for eQTL Browser
#
# Nick Burns
# June 2016

database <- function (conn = NULL) {
    
    init_ <- function () {
        drv <- DBI::dbDriver('MySQL')
        conn <- RMySQL::dbConnect(drv, dbname = 'eQTL_dw', default.file = "~/Documents/guest_db.cnf")
            
        return (conn)
    }

    close_ <- function() {
        RMySQL::dbDisconnect(conn)
    }

    if (missing(conn)) init_()
    else close_()
}

# parse_snp: given an RSID to query, return the chromosome and positions
# NOTE: that this uses the GLIDA package
parse_snp <- function (rsid) {
    return (glida::queryUCSC(glida::updatePositions(rsid)))
}

# browse: Queries the eQTL_dw
#
# Parameters:
#    target: string <- input$txt_query
#            a gene name or rsid 
#    type: string <- input$radio_query
#          specifies whether gene or rsid passed in
# Output:
#    data.table of resulting data
browse <- function (target, type = "gene") {
    
    query <- function () {
        
        lcl_query <- if (type == "gene") {
                sprintf("
                        SELECT
                            dimGene.gene_symbol,
                            dimTissue.SMTS,
                            factQTL.chromosome,
                            factQTL.build_37_pos,
                            factQTL.A1,
                            factQTL.A2,
                            factQTL.pvalue
                        FROM factQTL 
                            INNER JOIN dimGene ON factQTL.gene_id = dimGene.gene_id
                            INNER JOIN dimTissue on factQTL.tissue_id = dimTissue.tissue_id
                        WHERE dimGene.gene_symbol = '%s';", target)
        } else {
            ## need to do more here to turn RSID into chr, pos
            snp_info <- parse_snp(target)
            lcl_chr <- gsub("chr", "", snp_info$CHR)
            lcl_pos <- snp_info$POS
            
            print(snp_info)
            sprintf("SELECT * 
                    FROM factQTL
                    WHERE chromosome = %s
                      AND build_37_pos = %s", lcl_chr, lcl_pos)
        }
        
        return (lcl_query)
    }
    
    conn <- database()
    results <- data.table(dbGetQuery(conn, query()))
    database(conn)
    
    return (results)
}

dummy_data <- function (gene, data_type = "genotype_distributions") {
    
    
    tmp <- if (data_type == "genotype_distributions") {
        data.frame(genotype = rep(1:3, each = 100),
                   expression = c(rnorm(100, mean = 10, sd = 3),
                                  rnorm(100, mean = 15, sd = 4),
                                  rnorm(100, mean = 20, sd = 2.5)),
                   gene = gene)
    } else {
        data.frame(pvalue = sample(rpois(1000, 10)), 
                   build_37_pos = 1:1000,
                   gene = gene)
    }
    
    return (tmp)
}
