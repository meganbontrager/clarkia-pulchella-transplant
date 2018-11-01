Code and data to accompany the manuscript "Gene flow improves fitness at a range edge under climate change".  

by Megan Bontrager* and Amy L. Angert

Corresponding author: mgbontrager@gmail.com

Available on bioRxiv: https://doi.org/10.1101/399469

Analyses were conducted in R, version 3.4.4.

DIRECTORY/FILE:					DESCRIPTION:									
analysis_scripts                        	code to run most analyses in the paper
    	climate_fitness_effects.R           	referenced in figure 3, table S3
    	compare_local_foreign.R             	referenced in supp analyses 1, figure S5
    	gene_flow_effects.R                 	referenced in figure 4, table S4
    	genetic_differentiation_effects.R   	referenced in figure 5, table S5, table S6
    	local_vs_foreign_with_climate.R     	referenced in supplementary analyses 2
    	local_vs_foreign.R                  	referenced in figure S4, table S2

data                                    	clean data file used for analyses
    	Bontrager_transplant_data.csv       	all data used in main analyses, metadata supplied in Bontrager_transplant_metadata
    	Bontrager_transplant_metadata       	metadata for Bontrager_transplant_data

results
    	published_tables                    	tables arranged for publication, derived from raw tables
        	Bontrager_tableS1.csv
        	Bontrager_tableS2.csv
        	Bontrager_tableS3.csv
        	Bontrager_tableS4.csv
        	Bontrager_tableS5.csv
        	Bontrager_tableS6.csv
        	Bontrager_tableS7.csv
    	raw_tables                          	raw output of analysis scripts
        	fst_table .csv                  from script "genetic_differentiation_effects.R"
        	gf_table.csv                    from script "gene_flow_effects.R "
       		la_table_climate.csv            from script "local_vs_foreign_with_climate.R"
        	la_table.csv                    from script "local_vs_foreign.R”
		wi_table.csv			from script “climate_fitness_effects.R”
        
        
        
        
        
