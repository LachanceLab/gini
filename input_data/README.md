# Input data
Directory containing files and tables needed for other code to run properly. Below is a list of files needed, their filenames, and where to get them if not included in the repository already.

| Filename | Description | Source |
|:--------:|:--------:|:--------:|
| `chr_max_bps.txt` | List of the maximum base pairs per chromosome for binning purposes | This github repository |
| `traits_list.txt` | List of every trait (and more) used in the analysis, containing the trait description and shortened label, as well as the respective panUKB GWAS summary file data pulled from the phenotype manifest file | This github repository |
| `phenotype-description.csv` | List of every trait used by Prive et al. with trait descriptions and Prive et al.'s trait groups | [Prive et al. Github](https://github.com/privefl/UKBB-PGS/blob/main/phenotype-description.xlsx) (must convert to .csv file) |
| `phenotype-info.csv` | List of every trait used by Prive et al. with case/control information and heritability calculation statistics| [Prive et al. Github](https://github.com/privefl/UKBB-PGS/blob/main/phenotype-info.csv) |
| `pred-cor-PLR.csv` | List of every trait used by Prive et al. with polygenic score accuracy (partial correlation; pcor) information for each population | [Prive et al. Figshare](https://figshare.com/articles/dataset/Effect_sizes_for_215_polygenic_scores/14074760/2?file=31619357) |
| `PGS-effects-PLR.csv` | Giant matrix containing the effect size of every significant SNP with every trait analyzed by Prive et al. | [Prive et al. Figshare](https://figshare.com/articles/dataset/Effect_sizes_for_215_polygenic_scores/14074760/2?file=31619351) (must unzip first)|
| aau1043_datas3 | Recombination rate genome map created by Halldorsson et al. | Supplemental file of [Halldorsson et al. 2019](https://doi.org/10.1126/science.aau1043) (must be unzipped)