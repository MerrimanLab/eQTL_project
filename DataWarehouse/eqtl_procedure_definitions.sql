use eQTL_dw;

# reset staging table
drop procedure if exists qtl_reset_staging;
delimiter //
create procedure qtl_reset_staging()
begin
	truncate eqtl_staging;
end //
delimiter ;

# populate eQTL fact table
drop procedure if exists qtl_populate_fact_qtl;
delimiter //
create procedure qtl_populate_fact_qtl(IN lcl_tissue tinyint, IN lcl_source tinyint)
begin

	insert into factQTL (gene_id, tissue_id, chromosome, build_37_pos, A1, A2, beta, tstat, pvalue, source_id)
		SELECT
			gene.gene_id,
            lcl_tissue as 'tissue_id',
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
                pvalue,
                lcl_tissue as 'tissue_id'
			from eqtl_staging
            where pvalue < 0.01
        ) as stage
        inner join (
			select distinct ensembl_id, gene_id
            from dimGene
        ) as gene on gene.ensembl_id = stage.ensembl_id;
end //
delimiter ;


# populate expression fact table
# join to dimGene to get gene_id
# include tissue_id, as passed in as function 
# include source_id, so we know where this data came from
drop procedure if exists qtl_populate_fact_expr;
delimiter //
create procedure qtl_populate_fact_expr(IN lcl_tissue tinyint, IN lcl_source tinyint)
begin

	insert into factExpression (gene_id, tissue_id, source_id, rpkm)
	select g.gene_id,
		   lcl_tissue as 'tissue_id',
           lcl_source as 'source_id',
           s.rpkm
    from (
			select 
				substring_index(ensembl_id, '.', 1) as ensembl_id,
                gene_symbol,
                rpkm
            from expr_staging
    ) as s
    inner join dimGene as g on s.gene_symbol = g.gene_symbol AND s.ensembl_id = g.ensembl_id;
	
end //
delimiter ;
