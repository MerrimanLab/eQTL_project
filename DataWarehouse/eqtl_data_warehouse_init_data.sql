use eQTL_dw;


# populate meta-stagign
LOAD DATA LOCAL INFILE 'gtex_metadata.csv' INTO TABLE meta_staging FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 LINES;

select distinct SMTS, SMTSD from meta_staging;

# populate and update dimTissue
INSERT INTO dimTissue (smts, smtsd)
        SELECT DISTINCT SMTS, SMTSD from meta_staging;
        
UPDATE dimTissue
	set smts = 'Stomach'
    where SMTSD = 'Stomach';

UPDATE dimTissue
	set smts = 'Esophagus'
    where SMTSD = 'Esophagus - Mucosa';

UPDATE dimTissue
	set smts = 'Skin'
    where SMTSD = 'Skin - Sun Exposed (Lower leg)';
    
SELECT * FROM dimTissue;

## populate dimGene
LOAD DATA LOCAL INFILE '/mnt/DataDrive/gEXPR_eQTL_Datasets/GTEXData/GeneExpression/gtex_genes.csv'
INTO TABLE gene_staging
FIELDS TERMINATED BY ',' ENCLOSED BY '"'
LINES TERMINATED BY '\n'
IGNORE 1 LINES;

SELECT * FROM gene_staging limit 5;

insert into dimGene (ensembl_id, gene_symbol, chromosome, start_pos, end_pos, gene_biotype) 
	select * from gene_staging;