## Code (Part 1): Generating the traits table

The code in this folder is used to create a master traits table containing all the information known about each trait for this paper. Part 2 of our code then runs analyses and visualizations directly from this table.

### Scripts and the order to run them in

1. `assign_ancestries_to_IIDs.R`: Assigns an ancestry to each individual in the sample using Prive et al's code.
2. `download_clump_panUKB.sh`: Bash script that extracts panUKB GWAS summary statistics for each trait, uses PLINK to LD prune for independent GWAS hits, and saves the file using the Rscript in `filter_GWAS_independent.R`.
3. `liftover_hg38_hg19.py`: Converts the genetic map in `../input_data/aau1043_datas3` from GRCh38 to GRCh37 to match existing data. Used in supplementary analyses
4. `create_table.R`: Creates a table containing all relevant trait information: PGS accuracy, heritability, Gini coefficient, recombination rate, portability index, and PGS divergence.
5. `encode_sampled_genotypes.sh`: Runs a PLINK command that encodes the genotypes of the sampled indviduals for the top 100 bin SNPs so that `calculate_divergence.R` can easily run
6. `calculate_divergence.R`: Calculates the PRS for each sampled individual and calculates the level of divergence in PRS using the ANOVA F-stat
7. `reformat_ldscores.sh`: Cleans up the raw full LD scores from PanUKBB found in `~/input_data/ldscores/`.
8. `calculate_traitLD.R`: calculates a weighted average LD score for SNPs associated with a trait, as well as metrics of population-differences in LD scores

Helper scripts (no need to run independently; inside `helper_functions/`):
- `helper_functions.R`: contains functions used by multiple scripts related to binning, heritability calculations, and gini computation
- `filter_GWAS_independent.R`: Rscript that filters panUKB GWAS data to just LD independent SNPs. Called by `download_clump_panUKB.sh`
- `winners_curse_functions.R`: contains functions used by `create_table.R` to correct GWAS effect sizes for the winners curse. Adapted from this [repository](https://github.com/cpalmer718/gwas-winners-curse).