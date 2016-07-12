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
browse_qtls <- function (target, dimTissue, type = "gene") {
    
    query <- function () {
        
        lcl_query <- if (type == "gene") {
                sprintf("
                        SELECT
                            dimGene.gene_symbol,
                            factQTL.chromosome,
                            factQTL.build_37_pos,
                            factQTL.A1,
                            factQTL.A2,
                            factQTL.pvalue,
                            factQTL.tissue_id
                        FROM factQTL 
                            INNER JOIN dimGene ON factQTL.gene_id = dimGene.gene_id
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
    
    setkey(results, "tissue_id")
    results <- results[dimTissue]
    results <- results[!is.na(gene_symbol)]
    
    results[, tissue_score := min(pvalue), tissue_id]
    tissue_rank <- results[, .(rank_score = min(tissue_score)), tissue_id][order(rank_score)]
    tissue_rank[, rank := 1:.N]
    
    setkey(tissue_rank, tissue_id)
    results <- results[tissue_rank]
    
    results[rank > 4, smts := "Other"]
    
    return (results[, .(gene_symbol, chromosome, build_37_pos, A1, A2, pvalue, smts, SMTSD)])
}

lookup_tissues <- function () {
    
    conn <- database()
    tissue_info <- data.table(dbGetQuery(conn, "SELECT * FROM dimTissue;"))
    database(conn)
    
    setkey(tissue_info, "tissue_id")
    
    return (tissue_info)
}

# NOTE: something really weird happended trying to join factExpression with dimTissue
#       the cardinality blew out enormously. Resorting to a merge in R here to work around this.
browse_expression <- function (target) {
    
    query <- function () {
        lcl_query <- sprintf("
                                SELECT 
                                    g.gene_symbol, 
                                    g.chromosome, 
                                    g.gene_biotype,
                                    f.rpkm,
                                    f.tissue_id
                                FROM factExpression as f
                                  INNER JOIN dimGene as g ON g.gene_id = f.gene_id
                                WHERE g.gene_symbol = '%s';
                             ", target)
        return (lcl_query)
    }
    conn <- database()
    results <- data.table(dbGetQuery(conn, query()))
    dimTissue <- data.table(dbGetQuery(conn, "select * from dimTissue;"))
    database(conn)
    
    setkey(results, tissue_id)
    setkey(dimTissue, tissue_id)
    
    results <- results[dimTissue]
    
    return (results[!is.na(gene_symbol)])
}

# qtl_network
# this will change slightly, but the idea is there.
# having played with this, it doesn't make sense to do this
# between eQTL loci, but instead it should take a ref locus
# of GWAS hits and display the eQTLs to all genes.
# Code to visualise this is in eQTL_Arcdiagrams.Rmd
qtl_network <- function (ref_locus) {
    query <- function () {
        lcl_query <- sprintf("
        SELECT 
            g1.gene_symbol, 
            (g1.start_pos + g1.end_pos) / 2 as 'gene_midpoint', 
            f1.chromosome, 
            f1.build_37_pos, 
            f1.pvalue
        FROM factQTL_new f1
          INNER JOIN dimGene g1 ON g1.gene_id = f1.gene_id
        WHERE f1.chromosome = %s
          AND g1.gene_symbol != '%s'
          AND f1.build_37_pos IN (%s) ;
        ", ref_locus[1, chromosome], ref_locus[1, gene_symbol], paste0(ref_locus[pvalue < 0.0001, build_37_pos], collapse = ", "))
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

get_genes <- function (data_) {
    
    gene_ <- unique(data_[, gene_symbol])
    chr_ <- unique(data_[, chromosome])
    data_[, POS := build_37_pos]
    
    genes_in_region <- glida::queryUCSC(
        glida::fromUCSCEnsemblGenes(chromosome = chr_,
                                    start = data_[, min(build_37_pos)],
                                    end = data_[, max(build_37_pos)])
    )
    genes_in_region <- genes_in_region[genes_in_region$geneType == "protein_coding", ]
    
    return (genes_in_region)
}
display_eqtls <- function (data_) {
    
    gene_ <- unique(data_[, gene_symbol])
    chr_ <- unique(data_[, chromosome])
    
    # Formatting and variable creation
    # These are niceties to simplify the plotting and the interactive on_click
    data_[, POS := build_37_pos]
    data_[, position := build_37_pos / 1000000]
    data_[, association := -log10(pvalue + 1e-20)]
    
    genes_in_region <- get_genes(data_)
    
    viz <- ggplot(data_, aes(x = position, y = association)) +
        geom_point(aes(shape = factor(smts), 
                       size = 1., 
                       alpha = sqrt(1 / (pvalue + 1e-20))),
                   colour = "dodgerblue") +
        ylab("-log10( pvalue )") + xlab(sprintf("CHR%s position (MB)", chr_)) + ggtitle(sprintf("%s locus", gene_)) +
        guides(size = "none", alpha = "none") +
        scale_shape_discrete(name = "Top-ranked tissues") +
        theme_minimal()
    
    viz <- glida::geneAnnotation(viz, genes_in_region)
    
    return (viz)
    
}

display_expression <- function (gene_) {
    
    ggplot(browse_expression(gene_), aes(x = SMTSD, y = rpkm, group = SMTSD)) +
        geom_boxplot(aes(colour = smts, fill = smts), alpha = 0.5) +
        theme_minimal() +
        xlab("") +
        ggtitle(sprintf("Gene Expression: %s", gene_)) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
    
}

# browse_by_snp
# given a SNP (CHr, Pos, SNP, strand), find all nearby genes wtih eQTLs.
browse_by_snp <- function (snp) {
    conn <- database()
    results <- data.table(dbGetQuery(conn, sprintf("
                                     select distinct g.gene_symbol
                                    from factQTL f
                                      inner join dimGene g on g.gene_id = f.gene_id
                                    where f.build_37_pos BETWEEN %s ANd %s
                                      and f.chromosome = %s;
                                     ", snp$POS - 50, snp$POS + 50, 
                                                   gsub("chr", "", snp$CHR))))
    database(conn)
    
    return (results)
}

all_snp_info <- function (snp) {
    conn <- database()
    results <- data.table(dbGetQuery(conn, sprintf("
                                     select f.chromosome, f.build_37_pos, g.gene_symbol, f.pvalue
                                    from factQTL f
                                      inner join dimGene g on g.gene_id = f.gene_id
                                    where f.build_37_pos = %s
                                      and f.chromosome = %s;
                                     ", snp$build_37_pos, snp$chromosome)))
    database(conn)
    
    results <- results[, .(p.value = min(pvalue)), by = c("chromosome", "build_37_pos", "gene_symbol")][order(p.value, decreasing = FALSE)]
    return (results)
}
