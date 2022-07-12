## Code (Part 1): Generating the traits table

The code in this folder is used to create a master traits table containing all the information known about each trait for this paper. Part 2 of our code then runs analyses and visualizations directly from this table.

### Scripts and the order to run them in

1. `assign_ancestries_to_IIDs.R`: Assigns an ancestry to each individual in the sample using Prive et al's code.
2. `create_ancestry_AFs.sh`: PLINK script that calculates allele frequencies per ancestry
3. `append_betas_AFs.R`: Appends allele frequencies onto Prive et al.'s effect weights and generates files
4. `liftover.py`: Converts the genetic map in `../other_data/aau1043_datas3` from GRCh38 to GRCh37 to match existing data
5. `create_table.R`: Creates a table containing all relevant trait information, such as prediction performance, heritability, gini, portability, and PRS divergence.
6. `encode_sampled_genotypes.sh`: Runs a PLINK command that encodes the genotypes of the sampled indviduals for the top 100 bin SNPs so that `calculate_divergence.R` can easily run
7. `calculate_divergence.R`: Calculates the PRS for each sampled individual and calculates the level of divergence in PRS using the ANOVA F-stat