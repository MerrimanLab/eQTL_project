# Schema definition for eQTL data warehouse
#
# Nick Burns
# created: June 2016

drop database if exists eQTL_dw;
create database if not exists eQTL_dw;
use eQTL_dw;

# Staging tables
# --------------
# Used for bulk insert of raw GTEx data files.
create table if not exists eqtl_staging 
(
	ensembl_id varchar (32),
    snp_id varchar(32),
    beta float,
    tstat float,
    pvalue float
);
# Pvalues will be filtered to pvalue <= 0.01 (see data loading procedures for documentation)
create index idx_qtl_staging on eqtl_staging (pvalue);

# staging table to expression data
drop table if exists expr_staging;
create table if not exists expr_staging
(
	ensembl_id varchar(32),
    gene_symbol varchar(32),
    sample_id varchar(52),
    rpkm float
);
create index idx_expr_staging on expr_staging (ensembl_id, gene_symbol);

# staging table for metadata
drop table if exists meta_staging;
create table if not exists meta_staging
(
	sample_id varchar(32),
    smts varchar(36),
    smtsd varchar(36)
);

drop table if exists gene_staging;
create table gene_staging 
(
    ensembl_id varchar(32) not null,
    gene_symbol varchar(32) not null,
    chromosome tinyint not null,
    start_pos int not null,
    end_pos int not null,
    gene_biotype varchar(32)
);

# Data tables
# -----------
#
# Very simple schema:
#    - dimDataSource
#        * description of where the data came from (initially just GTEx, but will add to over time)
#    - dimGene
#        * all distinct ensembl_ids have been previously extracted
#        * ensembl_ids were mapped to gene symbols, chromosome, start and end positions using R (biomaRt)
#    - dimTissue
#        * eQTLs are generally tissue-specific, this maps each eQTL to the relevant tissue
#    - factQTL
#        * the main data table, contains QTLs and relevant statistics
# NOTE: the GTEx data is ~ 400 GB over 6 billion rows. To reduce this, we have filtered on pvalues and kept 
#       data types to the smallest possible types.

# dimGene
#    - contains enesmbl_id, gene symbol and genomic coordinates
#    - Genomic coordinates and biotype queried using biomaRt (see etl_gene_expression.Rmd)
#    - X, Y, Mt chromosomes ommited, as not analysed by GTEx with regards to eQTLs.
drop table if exists dimGene;
create table dimGene 
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
# the following index facilitates user queries
create index idx_gene_symbol on dimGene (gene_symbol, gene_id) using btree;
# the following index is inplace for the bulk load of factQTL,
# which is done from ensembl_id
create index idx_gene_ensembl on dimGene (ensembl_id, gene_id) using hash;
# index for bulk load of expression table
create index idx_gene_ids on dimGene (ensembl_id, gene_symbol, gene_id) using hash;

# dimTissue
#    - eQTLs are tissue-specific, this table facilitates that mapping
drop table if exists dimTissue;
create table if not exists dimTissue
(
	tissue_id tinyint not null auto_increment,
    smts varchar(52),
    smtsd varchar(64),
    primary key (tissue_id)
);

create table dimDataSource
(
	source_id tinyint not null auto_increment,
    source_name varchar(32) not null,
    source_description varchar(256),
    primary key (source_id)
);
insert into dimDataSource values (DEFAULT, "GTEx", "GTEx public release v6");

# factQTL
#    - main data table with all relevant eQTL statistics
#    - pvalues filtered to only keep pvalues <= 0.001 (see procedure definitions).
#    - alleles were initially omitted, but included now to facilitate joint GWAS:eQTL analysis
#        proposed by Zhu et al (Nature, 2016) 
#        (Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets)
#    - NOTE: that we cannot have a FK references to dimGene using gene_id - referential integrity needs 
#            to be imposed in the data load code!
create table factQTL 
(
	gene_id int not null,
    tissue_id tinyint not null,
    chromosome tinyint not null,
    build_37_pos int not null,
    A1 varchar(24),
    A2 varchar(24),
    beta float,
    tstat float,
    pvalue float,
    source_id tinyint not null,
    foreign key (tissue_id) references dimTissue (tissue_id),
    foreign key (source_id) references dimDataSource (source_id)
);
#create index idx_qtl_gene_tissue on factQTL (gene_id, tissue_id);  # query based on gene-expression
# create index idx_qtl_snp on factQTL (chromosome, build_37_pos);   # query based on snp

# factExpression  
# Gene Expression data  
drop table if exists factExpression;
create table factExpression 
(
	gene_id int not null,
    tissue_id tinyint not null,
    source_id tinyint not null,
    rpkm float not null,
    foreign key (gene_id) references dimGene (gene_id),
    foreign key (tissue_id) references dimTissue (tissue_id),
    foreign key (source_id) references dimDataSource (source_id)
);
create index idx_expression on factExpression (gene_id, tissue_id) using hash;