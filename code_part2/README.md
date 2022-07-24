## Code (Part 2): Generating the traits table

The code in this folder is used to analyze and visualize the data generated in Part 1 of the code, particularly the traits table. Also contains supplementary analyses.

### Scripts (no order required to run scripts unless specified)

- `plot_divergence.R`: plots the PRS distribution of different ancestries, highlighting their divergence
- `plot_lorenz_curve.R`: plots a Lorenz curve for a trait using its significant bins' genetic variance contribution
- `plot_scatterplot_matrix_density.R`: plots a scatterplot matrix of the main six statistics and the densities of those statistics
- `plot_statistic_PCA.R`: plots a PCA plot of the 211 traits, using the main six summary statistics as the initial dimensions
- `generate_extra_ginis.R`: generates Gini coefficients using different calculation settings and then analyzes the difference in settings

Helper scripts (no need to run independently):
- `custom_ggbiplot.R`: modified version of the `ggbiplot` function of the `ggbiplot` package that allows customization of certain aesthetic plot attributes