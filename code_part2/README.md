## Code (Part 2): Generating the traits table

The code in this folder is used to analyze and visualize the data generated in Part 1 of the code, particularly the traits table. Also contains supplementary analyses.

### Scripts (no order required to run scripts unless specified)

- `plot_gini_m_D.R`: plots Figure 2 of paper (high and low Gini, portability, and PGS divergence examples).
- `plot_scatterplot_matrix_density.R`: plots a scatterplot matrix of the main six statistics, the densities of those statistics, and the scatterplot between prevalence and the six statistics. Figures 3-5 of paper.
- `plot_statistic_PCA.R`: plots a PCA plot of all traits, using the main six summary statistics as the initial dimensions. Plots both main paper Figure 5 and supplemental figure S8.
- `gini_PLR_vs_LDpred2.R`: compares Gini coefficients when using effect sizes derived from PLR versus LDpred2. Figure S4. Requires downloading new Prive data.
- `gini_supplemental.R`: generates Gini coefficients using different calculation settings and then analyzes the difference in settings, showing robustness of Gini. Supplemental figures S5-S7.
- `print_tables.R`: generates and prints tables in main paper (Table 1 and 2) and supplemental (Table S1)

Helper scripts (no need to run independently):
- `custom_ggbiplot.R`: modified version of the `ggbiplot` function of the `ggbiplot` package that allows customization of certain aesthetic plot attributes for our paper.