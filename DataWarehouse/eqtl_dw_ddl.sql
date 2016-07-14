# eQTL Datawarehouse
# Data definition language (DDL) script
#
# Contains the specific DDL commands to establish the eQTL DW.
# The schema (and therefore, this script) is organised as follows:
# 
#     Section 1: Staging
#     ------------------
#         - contains relevant staging tables for bulk data loads and 
#           intermediate transformations
#         - tables: stage1_qtl, stage2_qtl, stage_expr, stage_meta, stage_gene
#
#     Section 2: Dimension tables  
#     ---------------------------
#         - tables: dim_gene, dim_tissue, dim_data_source
#         - indexes also defined on these tables
#
#     Section 3: Fact tables
#     ----------------------
#         - tables: 
#              fact_qtl: main data table with relevant eQTL statistics,
#                        pvalues filtered (p < 0.1) prior to bulk load
#              fact_expr: main data table with gene expression statistics
#         - ideally, these two fact tables would be one and the same 
#           (fact_genome_stats, or something similar). However, for optimal efficiency
#           and also convenience, we have separated these out for this application.
#         - relevant indexes also defined
#
#     Section 4: Stored procedures
#     ----------------------------
#         - ETL procedures: the majority of the ETL processes are encoded in procs
#         - APP procedures: as the app matures, the main queries will be encoded in procs.
#
# Nick Burns
# July 2016


# 
# Database creation
#
drop database if exists qtldw;
create database qtldw;
use qtldw;

#
# Section 1: Staging
#

# stage1_qtl:  initial landing table for bulk load of GTEx QTL files
# ----------
drop table if exists stage1_qtl;
create table if not exists stage1_qtl 
(
    ensembl_id varchar (32),
    snp_id varchar(128),
    beta float,
    tstat float,
    pvalue float
);

# stage2_qtl: intermediate table between stage1_qtl and fact_qtl
# -----------
#    - the whole purpose of this table is to reduce the impact of the string parsing on stage1_qtl.snp_id
#    - A1 and A2 are necessarily large, as there are large INDELs in the GTEx data.
#    - Index on ensembl_id to facilitate join between stage2_qtl and dim_gene
#      during unstaging.
drop table if exists stage2_qtl;
create table if not exists stage2_qtl
(
    ensembl_id varchar(32),
    chromosome tinyint,
    build_37_pos int,
    A1 varchar(128),
    A2 varchar(128),
    beta float,
    tstat float,
    pvalue float,
    tissue_id int,
    source_id int
);
create index idx_unstage_qtl on stage2_qtl (ensembl_id);

# stage_expr:  staging table for bulk load of GTEx gene expression files
drop table if exists stage_expr;
create table if not exists stage_expr
(
    ensembl_id varchar(32),
    gene_symbol varchar(32),
    sample_id varchar(52),
    rpkm float 
);
create index idx_stage_expr on stage_expr (ensembl_id, gene_symbol);

# stage_meta: staging table for bulk load of GTEx sample / tissue metadata
drop table if exists stage_meta;
create table if not exists stage_meta
(
    sample_id varchar(32),
    smts varchar(128),
    smtsd varchar(128)
);

# stage_gene: staging table for bulk load of GTEx gene information
drop table if exists stage_gene;
create table if not exists stage_gene
(
    ensembl_id varchar(32),
    gene_symbol varchar(32),
    chromosome tinyint,
    start_pos int,
    end_pos int,
    gene_biotype varchar(32)
);

#
# Section 2: Dimension tables
#

# dim_gene:
# ---------
#   Gene information taken originally from GTEx gene expression data
#   Biotype and genomic coordinates added using biomaRt (see etl_gene_expression.Rmd)
#   X, Y, Mt chromosomes omitted, as not analysed by GTEx wtih regards to eQTLs.
drop table if exists dim_gene;
create table dim_gene 
(
    gene_id int auto_increment not null,
    ensembl_id varchar(32) not null,
    gene_symbol varchar(32) not null,
    chromosome tinyint not null,
    start_pos int not null,
    end_pos int not null,
    gene_biotype varchar(32),
    primary key (gene_id)
);
# idx_ensembl: this index is purely for population of fact_qtl which joins to ensembl_ids
# idx_gene_ids: this index is purely for population of fact_expr
# idx_gene_symbol facilitate App lookups.
create index idx_ensembl on dim_gene (ensembl_id, gene_id) using hash;
create index idx_gene_ids on dim_gene (ensembl_id, gene_symbol, gene_id) using hash;
create index idx_gene_symbol on dim_gene (gene_symbol, gene_id) using btree;

# dim_tissue:
# -----------
#   Tissue information from GTEx.
#     - smts is a 'general categorisation' of tissues e.g. brain, adipose, blood etc...
#     - smtsd is a more specific categorisation e.g. 'Brain cerebellum', 'Brain cortex', etc.
#   This table will never exceed a few dozen rows (currently ~ 50 rows), so no need for indexing at this stage.
drop table if exists dim_tissue;
create table if not exists dim_tissue
(
    tissue_id tinyint not null auto_increment,
    smts varchar(128),
    smtsd varchar(128),
    primary key (tissue_id)
);

# dim_data_source: lookup table to record where data was originally extracted from
drop table if exists dim_data_source;
create table if not exists dim_data_source
(
    source_id tinyint not null auto_increment,
    source_name varchar(32) not null,
    source_description varchar(256),
    primary key (source_id)
);
insert into dim_data_source values (DEFAULT, "GTEx", "GTEx public release v6 2016");

#
# Section 3: Fact tables
#

# fact_qtl: main data table with all relevant eQTL statistics
drop table if exists fact_qtl;
create table if not exists fact_qtl 
(
    gene_id int not null,
    tissue_id tinyint not null,
    chromosome tinyint not null,
    build_37_pos int not null,
    A1 varchar(128),
    A2 varchar(128),
    beta float,
    tstat float,
    pvalue float,
    source_id tinyint not null,
    foreign key (tissue_id) references dim_tissue (tissue_id),
    foreign key (source_id) references dim_data_source (source_id)
);
# idx_by_gene: index to facilitate querying by gene (including genomic coordinates)
# idx_by_coord: index to facilitate query by coordinate (including gene)
create index idx_by_gene on fact_qtl (gene_id, chromosome, build_37_pos) using hash;
create index idx_by_coord on fact_qtl (chromosome, build_37_pos, gene_id) using btree;

# fact_expr: main data table with GTEx gene expression data   
drop table if exists fact_expr;
create table fact_expr 
(
    gene_id int not null,
    tissue_id tinyint not null,
    source_id tinyint not null,
    rpkm float not null,
    foreign key (gene_id) references dim_gene (gene_id),
    foreign key (tissue_id) references dim_tissue (tissue_id),
    foreign key (source_id) references dim_data_source (source_id)
);
# idx_expression: query expression by gene, and optionally by tissue
create index idx_expression on fact_expr (gene_id, tissue_id) using hash;

#
# Section 4: Stored procedures
# 

# udf_parse_qtl:  parse qtls from stage1_qtl -> stage2_qtl
# Note: have placed the parameters lcl_tissue and lcl_soruce here
#          to give udf_unstage_qtl the best chance of maximising indexes
#          In SQL Server, parametised functions like these cannot be
#          cached, and I just dont know how MySQL / MariaDB handles this.
drop procedure if exists udf_parse_qtl;
delimiter //
create procedure udf_parse_qtl(IN lcl_tissue tinyint, IN lcl_source tinyint)
begin
    insert into stage2_qtl
    select 
        substring_index(ensembl_id, '.', 1) as ensembl_id,
        substring_index(snp_id, '_', 1) as chromosome,
        substring_index(substring_index(snp_id, '_', -4), '_', 1) as build_37_pos,
        substring_index(substring_index(snp_id, '_', -3), '_', 1) as A1,
        substring_index(substring_index(snp_id, '_', -2), '_', 1) as A2,
        beta,
        tstat,
        pvalue,
        lcl_tissue as 'tissue_id',
        lcl_source as 'source_id'
    from stage1_qtl;
end //
delimiter ;

# udf_unstage_qtl: unstage qtls from stage2_qtl -> fact_qtl
drop procedure if exists udf_unstage_qtl;
delimiter //
create procedure udf_unstage_qtl()
begin
    insert into fact_qtl (gene_id, chromosome, build_37_pos, A1, A2, beta, tstat, pvalue, tissue_id, source_id)
    select 
        g.gene_id,
        q.chromosome,
        q.build_37_pos,
        q.A1,
        q.A2,
        q.beta,
        q.tstat,
        q.pvalue,
        q.tissue_id,
        q.source_id
    from stage2_qtl as q
      inner join dim_gene as g on g.ensembl_id = q.ensembl_id;
end //
delimiter ;


# udf_unstage_expr: unstage expression data from stage_expr -> fact_expr
drop procedure if exists udf_unstage_expr;
delimiter //
create procedure udf_unstage_expr(IN lcl_tissue tinyint, IN lcl_source tinyint)
begin

    insert into fact_expression (gene_id, tissue_id, source_id, rpkm)
    select 
        g.gene_id,
	lcl_tissue as 'tissue_id',
        lcl_source as 'source_id',
        s.rpkm
    from (
	select 
	    substring_index(ensembl_id, '.', 1) as ensembl_id,
            gene_symbol,
            rpkm
        from stage_expr
    ) as s
    inner join dim_gene as g 
      on s.gene_symbol = g.gene_symbol AND s.ensembl_id = g.ensembl_id;	
end //
delimiter ;
