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

## Create INDEX {gene, chromosome, tissue}
## Let most likely query will be: "for a given gene, which is on chromosome C,
## return all eQTLs (filtered / grouped by tissue)
create index idx_eqtls on factQTL (ensembl_id, tissue) using btree;
## note {ensembl_id, chromosome} is unqiue anyway, so no need to specify chromosome.
## eQTLs are all defined relative to a GENE, so a CHR:start-end lookup is not relevant.


drop procedure if exists populateFact;

delimiter //
create procedure populateFact(IN lcl_tissue int)
begin
	   
	INSERT INTO factQTL (ensembl_id, tissue, chromosome, build_37_pos, beta, tstat, pvalue, source_name)
		SELECT 
			substring_index(ensembl_id, '.', 1),
			lcl_tissue,
			substring_index(snp_id, '_', 1) as chromosome,
			substring_index(substring_index(snp_id, '_', 2), '_', -1),
			beta,
			tstat,
			pvalue,
			'GTEx'
		FROM eQTL_staging ;
    
END //
delimiter ;
