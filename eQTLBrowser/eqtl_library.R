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
# NOTE: assuming GWAS colnames (CHR, POS, P)
qtl_network <- function (chr_, start_, end_) {
    query <- function () {
        lcl_query <- sprintf("
        SELECT 
            g1.gene_symbol, 
            (g1.start_pos + g1.end_pos) / 2 as 'gene_midpoint', 
            f1.chromosome, 
            f1.build_37_pos, 
            f1.pvalue
        FROM factQTL f1
          INNER JOIN dimGene g1 ON g1.gene_id = f1.gene_id
        WHERE f1.chromosome = %s
          AND f1.build_37_pos BETWEEN %s AND %s ;
        ", chr_, start_, end_)
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
    
    chr_ <- unique(data_[, chromosome])
    
    if (("build_37_pos" %in% colnames(data)) & (! "POS" %in% colnames(data))) {
        data_[, POS := build_37_pos]
    }
    
    genes_in_region <- glida::queryUCSC(
        glida::fromUCSCEnsemblGenes(chromosome = chr_,
                                    start = data_[, min(POS)],
                                    end = data_[, max(POS)])
    )
    genes_in_region <- genes_in_region[genes_in_region$geneType == "protein_coding", ]
    
    return (genes_in_region)
}
display_eqtls <- function (data_, gwas_data = NULL,
                           show_genes = TRUE, show_tissues = TRUE, alpha_pvalues = TRUE,
                           show_title = TRUE) {
    
    gene_ <- unique(data_[, gene_symbol])
    chr_ <- unique(data_[, chromosome])
    
    # Formatting and variable creation
    # These are niceties to simplify the plotting and the interactive on_click
    data_[, POS := build_37_pos]
    data_[, position := build_37_pos / 1000000]
    data_[, association := -log10(pvalue + 1e-20)]
    
    # weight the size of points
    if (!is.null(gwas_data)) {
        data_[, gwas_weight_ := weight_(gwas_data, build_37_pos), by = build_37_pos]
    } else {
        data_[, gwas_weight_ := 1]
    }
    
    
    viz <- ggplot(data_, aes(x = position, y = association)) +
        geom_point(aes(shape = ifelse(show_tissues, factor(smts), 'a'),
                       size = gwas_weight_),
                   alpha = ifelse(alpha_pvalues, sqrt(1 / (data_[,pvalue] + 1e-20)), 0.3),
                   colour = "darkblue") +
        ylab("-log10( pvalue )") + xlab(sprintf("CHR%s position (MB)", chr_)) + 
        ggtitle(ifelse(show_title, sprintf("%s locus", gene_), "")) +
        guides(size = "none", alpha = "none", shape = ifelse(show_tissues, NULL, "none")) +
        theme_minimal()
    
    if (show_genes) {
        genes_in_region <- get_genes(data_)
        viz <- glida::geneAnnotation(viz, genes_in_region)
    }
    if (show_tissues) {
        viz <- viz + scale_shape_discrete(name = "Top-ranked tissues")
    }
    
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
                                    where f.build_37_pos BETWEEN %s AND %s
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


#### GWAS : QTL Functions  

extract_gwas <- function (file_, chr_, start_, end_) {
    
    gwas_ <- fread(file_)
    
    return (gwas_[(CHR == chr_ & POS >= start_ & POS <= end_)])
}

display_gwas <- function (data_) {
    
    chr_ <- unique(data_[, CHR])
    data_[, chromosome := CHR]
    
    # Formatting and variable creation
    # These are niceties to simplify the plotting and the interactive on_click
    data_[, position := POS / 1000000]
    data_[, association := -log10(P + 1e-20)]
    
    genes_in_region <- get_genes(data_)
    
    viz <- ggplot(data_, aes(x = position, y = association)) +
        geom_point(colour = "dodgerblue", alpha = 0.3) +
        ylab("-log10( pvalue )") + xlab("position (MB)") + 
        ggtitle(sprintf("Chromosome %s : %s MB - %s MB", chr_, min(data_[, position]), max(data_[, position]))) +
        theme_minimal()
    
    return (viz)
    
}

display_qtl_network <- function (viz, long_range_qtls, show_endpoint = TRUE) {
    print("trying to display the qtl network...")
    layer_ <- viz +
        geom_curve(data = long_range_qtls,
                   aes(gene_midpoint / 1000000, xend = build_37_pos / 1000000, 
                       y = -10, yend = -log10(pvalue + 1e-20),
                       alpha = 1 / (pvalue + 1e-20)),
                   colour = "darkgrey", curvature = 0.3) +
        geom_text(data = unique(long_range_qtls[, .(gene_symbol, gene_midpoint)]),
                  aes(x = gene_midpoint / 1000000, y = -10, label = gene_symbol)) +
        guides(alpha = "none") +
        theme_minimal()
    
    return (layer_)
}

nearest_neighbours <- function (ref_position, gwas) {
    
    dx <- as.matrix(dist(c(ref_position, gwas$POS)))
    dx <- sort(dx[1, ], decreasing = FALSE)
    
    return (dx)
}
weighted_pvalue <- function (gwas, dx) {
    sum(-log10(gwas$P + 1e-50) * (1 / dx))
}

weight_ <- function (gwas, x) {
    
    nearest_ <- function () {
        DX <- as.matrix(dist(c(x, gwas$POS)))[1, -1]
        
        return (DX[order(DX, decreasing = FALSE)])
    }
    size_ <- sum(-log10(gwas$P) / nearest_())
    
    return (size_)
}

