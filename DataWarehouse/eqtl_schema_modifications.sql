# Updates to the eQTL_dw schema
#
# Have populated the data warehouse. The data size is ~ 410 GB, and the index size is ~ 350 GB.
# I need to trim this back a little. The plan is to:
#
#    - recode dimGene (ensembl_id, gene_symbol, chromosome)
#         to: dimGene (gene_id, ensembl_id, gene_symbol, chromosome, start, end)
#
#    - change the foreign key from factQTL -> dimGene to be on the int (gene_id)
#    -
#
#
# references: datatypes storage size: http://dev.mysql.com/doc/refman/5.7/en/storage-requirements.html

use eQTL_dw;

alter table dimGene rename dimGene_old;
alter table dimTissue rename dimTissue_old;
alter table factQTL rename factQTL_old;

# NOTE: can't have a primary key here, because not all ensembl_ids mapped to a gene symbol :(

create table dimGene 
(
	gene_id smallint not null,
    ensembl_id varchar(32) not null,
    gene_symbol varchar(16),
    chromosome tinyint not null,
    start_pos int,
    end_pos int,
    primary key (ensembl_id, gene_symbol)
);
create index idx_gene on dimGene (gene_id);

create table dimTissue
(
	tissue_id tinyint not null auto_increment,
    tissue_description varchar(128) not null,
    primary key (tissue_id)
);
insert into dimTissue select * from dimTissue_old;

create table factQTL 
(
	gene_id smallint not null,
    tissue tinyint not null,
    chromosome tinyint not null,
    build_37_pos int not null,
    beta float,
    tstat float,
    pvalue float,
    source_name varchar(32) not null,
    foreign key (gene_id) references dimGene (gene_id),
    foreign key (tissue) references dimTissue (tissue_id),
    foreign key (source_name) references dimDataSource (source_name)
);
create index idx_qtl on factQTL (gene_id, tissue);


# The data is currently sitting in factQTL_old
# I need to get this into factQTL, but have to do it in chunks
# because it is too large to do in one go. Here we go.
drop procedure if exists populateFact;
delimiter //
create procedure populateFact()
begin
	   
	declare lcl_tissue int;
    declare chunk int;
    set lcl_tissue = 1;
    
    while (lcl_tissue < 3)
    do
		
		set chunk = 0;
        
        while (chunk < 10)
        do
			SELECT lcl_tissue as 'current tissue', chunk as 'current chunk', current_time() as 'time';
			INSERT INTO factQTL (gene_id, tissue, chromosome, build_37_pos, beta, tstat, pvalue, source_name)
				SELECT 
					g.gene_id, f.tissue, g.chromosome, f.build_37_pos, f.beta, f.tstat, f.pvalue, f.source_name
				FROM factQTL_old f
					INNER JOIN (
						SELECT distinct gene_id, ensembl_id, chromosome from dimGene
					) g ON g.ensembl_id = f.ensembl_id
				WHERE f.tissue = lcl_tissue 
                  AND g.chromosome % 10 = chunk
                  AND f.pvalue < 0.1;
                
			set chunk = chunk + 1;
		end while;
            
		set lcl_tissue = lcl_tissue + 1;
	end while ;   
END //
delimiter ;

drop procedure if exists filterPvalues;
delimiter //
create procedure filterPvalues()
begin
	   
	declare lcl_tissue int;
    declare chunk int;
    set lcl_tissue = 1;
    
    while (lcl_tissue < 44)
    do
		
		set chunk = 0;
        
        while (chunk < 10)
        do
			delete from factQTL_old where pvalue > 0.1;
                
			set chunk = chunk + 1;
		end while;
            
		set lcl_tissue = lcl_tissue + 1;
	end while ;   
END //
delimiter ;


show processlist;

