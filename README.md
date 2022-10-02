Data and analysis from: "Body mass, temperature, and depth shape productivity in sharks and rays"
===========
  
This repository contains the data and model-fitting process used for the article
"Body mass, temperature, and depth shape productivity in sharks and rays", which
is currently in press at *Ecology and Evolution*. 

## Description of data and file structure

The dataset includes life history information and phylogeny for 63 species of 
sharks(40 spp), rays (20 spp), and chimaeras (3 spp), and contains three files:

`rmax-scaling-analysis.R`: R code with minimum working example of how to load 
data files, fit models phylogenetic linear models using the `pgls` function in
the `caper` package, run information-theoretic comparisons, and check diagnostics.

`chond-data.csv`: Data frame with 63 rows (species) and 11 variables. Some of 
these variables are based on the same life history trait but are transformed for
ease of interpretation and analysis.

The 11 variables are:
- `id` (int): ordinal variable 
- `species`(chr): scientific name  
- `maxwt_all` (num): maximum weight, in grams
- `log_wt` (num):  natural log of maximum weight
- `rmax` (num): maximum intrinsic rate of population increase, *r<sub>max</sub>*
- `t_median` (num): median environmental temperature in Celsius
- `t_med_corr` (num):  median environmental temperature in Celsius, corrected for mesothermic species. This is the temperature used in analyses as it represents median metabolic temperature
- `invtemp` (num): inverse temperature, calculated as *1/k<sub>B</sub>T*, where *k<sub>B</sub>* is the Boltzmann constant (8.617 × 10 −5 eV) and *T* is median corrected temperature in Kelvin
- `invtemp_scaled` (num): scaled `invtemp` 
- `depth_median` (num): median depth, in metres
- `depth_scaled` (num): scaled `depth_median`

`stein-et-al-single.tree`: Phylogenetic tree with scaled branch lengths from Stein et al. (2018) used in analyses. These are freely downloadable from <http://vertlife.org/sharktree/>.

## Sharing/access Information

These data and analysis can also be found at <https://github.com/sebpardo/shark-rmax-scaling>.