# eQTL Datawarehouse
# Initial Data Population Script
#
# Loads lookup data into staging and then dimension tables. Specifically:
#     - file: gtex_metadata.csv -> stage_meta -> dim_tissue
#     - file: gtex_genes.csv -> stage_gene -> dim_gene
#
# NOTES:
#   1. This script should not be run automatically. Users should run each statement manually,
#      in particular, users should review the UPDATE dim_tissue statements to see if they need to
#      be run.
#   2. This script only needs to be run once. After running, back up the database.
#
# Nick Burns
# July 2016
use qtldw;


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
    
# Turns out that both the Stomach and esophagus entries are 
# double ups. Will remove both of these
delete from dimTissue where tissue_id IN (54, 55, 58);
    
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
