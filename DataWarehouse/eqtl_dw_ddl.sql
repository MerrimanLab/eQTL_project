create database if not exists eQTL_dw;
use eQTL_dw;

create table if not exists eQTL_staging 
(
	ensembl_id varchar (32),
    snp_id varchar(32),
    beta float,
    tstat float,
    pvalue float
);

create table if not exists dimGene
(
    ensembl_id varchar(32),
    gene_symbol varchar(32),
    chromosome varchar(2) not null,
    primary key (ensembl_id),
    constraint unique(ensembl_id, gene_symbol),
    constraint unique (ensembl_id, gene_symbol, chromosome)
);

create table if not exists dimTissue
(
	tissue_id int not null auto_increment,
    tissue_description varchar(128) not null,
    primary key (tissue_id)
);

create table if not exists dimDataSource
(
    source_name varchar(32) not null,
    source_description varchar(512) null,    # long-form description if required
    primary key (source_name)
);
insert into dimDataSource values ("GTEx", "GTEx public eQTL release V6.");

create table if not exists factQTL
(
	ensembl_id varchar(32) not null,
    tissue int not null,
    chromosome varchar(2) not null,
    build_37_pos int not null,
    beta float,
    tstat float,
    pvalue float,
    source_name varchar(32) not null,
    foreign key (ensembl_id) references dimGene (ensembl_id),
    foreign key (tissue) references dimTissue (tissue_id),
    foreign key (source_name) references dimDataSource (source_name)
);
