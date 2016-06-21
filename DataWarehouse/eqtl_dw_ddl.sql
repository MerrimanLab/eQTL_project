# Schema definition for eQTL data warehouse
#
# Nick Burns
# created: June 2016

create database if not exists eQTL_dw;
use eQTL_dw;

# Staging table
# -------------
# Used for bulk insert of raw GTEx data files.
create table if not exists eqtl_staging 
(
	ensembl_id varchar (32),
    snp_id varchar(32),
    beta float,
    tstat float,
    pvalue float
);
# Pvalues will be filtered to pvalue <= 0.001 (see data loading procedures for documentation)
create index idx_staging on eqtl_staging (pvalue);

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
#    - 9 ensembl_ids map to more than one gene symbol
#        e.g. ENSXXXX -> MIR50A, MIR50B, MIR50C
#        only seems to affect SNORA and MIR genes
#    - Given this restriction, the natural primary key is the
#        combination of (gene_symbol, ensembl_id).
#        However, users are most likely to query by gene_symbol and then
#        join to the fact table via gene_id. Thus, we will define a PK:
#        (gene_symbol, gene_id) to facilitate this lookup / join.
#        NOTE: this means all factQTL lookups should include:
#                 ... WHERE dimGene.gene_symbol = x
create table dimGene 
(
	gene_id smallint not null,
    ensembl_id varchar(32) not null,
    gene_symbol varchar(16) not null,
    chromosome tinyint not null,
    start_pos int not null,
    end_pos int not null
);
# the following index facilitates user queries
create index idx_gene_symbol on dimGene (gene_symbol, gene_id) using btree;
# the following index is inplace for the bulk load of factQTL,
# which is done from ensembl_id
create index idx_gene_ensembl on dimGene (ensembl_id, gene_id) using hash;

# dimTissue
#    - eQTLs are tissue-specific, this table facilitates that mapping
create table dimTissue
(
	tissue_id tinyint not null auto_increment,
    tissue_description varchar(128) not null,
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
	gene_id smallint not null,
    tissue_id tinyint not null,
    chromosome tinyint not null,
    build_37_pos int not null,
    A1 varchar(2),
    A2 varchar(2),
    beta float,
    tstat float,
    pvalue float,
    source_id tinyint not null,
    foreign key (tissue_id) references dimTissue (tissue_id),
    foreign key (source_id) references dimDataSource (source_id)
);
create index idx_qtl_gene_tissue on factQTL (gene_id, tissue_id);