## Code (Part 2): Generating the traits table

The code in this folder is used to analyze and visualize the data generated in Part 1 of the code, particularly the traits table. Also contains supplementary analyses.

### Scripts (no order required to run scripts unless specified)

- `plot_gini_m_D.R`: plots high and low example traits of select metrics.
- `plot_scatterplot_matrix_density.R`: plots a scatterplot matrix of the main six statistics and the densities of those statistics
- `plot_statistic_PCA.R`: plots a PCA plot of all traits, using the main six summary statistics as the initial dimensions.
- `gini_roustness.R`: generates Gini coefficients using different calculation settings and then analyzes the difference in settings, showing robustness of Gini.
- `print_tables.R`: generates and prints tables in main paper (Table 1 and 2) and supplemental (Table S1)

Helper scripts (no need to run independently):
- `custom_ggbiplot.R`: modified version of the `ggbiplot` function of the `ggbiplot` package that allows customization of certain aesthetic plot attributes for our paper.
