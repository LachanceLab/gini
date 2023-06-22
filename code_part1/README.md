## Code (Part 1): Generating the traits table

The code in this folder is used to create a master traits table containing all the information known about each trait for this paper. Part 2 of our code then runs analyses and visualizations directly from this table.

### Scripts and the order to run them in

0. `harmonize_1kG.sh`: Used to download 1000 Genomes Project Phase3 data, filtered to PanUKBB SNPs. Also generates list of samples by population as well as EUR-only samples.
1. `assign_ancestries_to_IIDs.R`: Assigns an ancestry to each individual in the sample using Prive et al's code.
2. `download_clump_panUKB.sh`: Bash script that extracts panUKB GWAS summary statistics for each trait, uses PLINK to LD prune for independent GWAS hits, and saves the file using the Rscript in `filter_GWAS_independent.R`.
3. `create_table.R`: Creates a table containing all relevant trait information: PGS accuracy, heritability, Gini coefficient, recombination rate, portability index, and PGS divergence.
4. `calculate_divergence.sh`: Computes PGS for all traits on a subset of the UKB population using PLINK2
5. `calculate_divergence.R`: Combines the chromosome-specific PGS data from `calculate_divergence.sh` and computes an ANOVA by ancestry on PGS for each trait
6. `reformat_ldscores.sh`: Cleans up the raw full LD scores from PanUKBB found in `~/input_data/ldscores/`.
7. `calculate_traitLD.R`: calculates a weighted average LD score for SNPs associated with a trait, as well as metrics of population-differences in LD scores

Helper scripts (no need to run independently; inside `helper_functions/`):
- `helper_functions.R`: contains functions used by multiple scripts related to gvc and gini computation
- `filter_GWAS_independent.R`: Rscript that filters panUKBB GWAS data to just LD independent SNPs. Called by `download_clump_panUKB.sh`
- `winners_curse_functions.R`: contains functions used by `create_table.R` to correct GWAS effect sizes for the winners curse. Adapted from this [repository](https://github.com/cpalmer718/gwas-winners-curse).