use eQTL_dw;

# reset staging table
drop procedure if exists qtl_reset_staging;
delimiter //
create procedure qtl_reset_staging()
begin
	delete from eqtl_staging;
end //
delimiter ;

# populate fact table
drop procedure if exists qtl_populate_fact;
delimiter //
create procedure qtl_populate_fact(IN lcl_tissue tinyint, IN lcl_source tinyint)
begin
	insert into factQTL (gene_id, tissue_id, chromosome, build_37_pos, A1, A2, beta, tstat, pvalue, source_id)
		SELECT
			gene.gene_id,
            lcl_tissue,
            stage.chromosome,
            stage.build_37_pos,
            stage.A1,
            stage.A2,
            stage.beta,
            stage.tstat,
            stage.pvalue,
            lcl_source		
        FROM (
			select 
				substring_index(ensembl_id, '.', 1) as ensembl_id,
                substring_index(snp_id, '_', 1) as chromosome,
                substring_index(substring_index(snp_id, '_', -4), '_', 1) as build_37_pos,
                substring_index(substring_index(snp_id, '_', -3), '_', 1) as A1,
                substring_index(substring_index(snp_id, '_', -2), '_', 1) as A2,
                beta,
                tstat,
                pvalue
			from eqtl_staging
            where pvalue < 0.001
        ) as stage
        inner join (
			select distinct ensembl_id, gene_id
            from dimGene
        ) as gene on gene.ensembl_id = stage.ensembl_id;
end //
delimiter ;




