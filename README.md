# MLLPT <img src='man/figures/logo.png' align="right" height="139" />
<!-- badges: start -->
[![Travis build status](https://travis-ci.com/T-Heide/MLLPT.svg?token=cqqKUEunNazVxSmjTbrG&branch=master)](https://travis-ci.com/T-Heide/MLLPT)
[![Codecov test coverage](https://codecov.io/gh/T-Heide/MLLPT/branch/master/graph/badge.svg?token=1M0X8H2FRY)](https://codecov.io/gh/T-Heide/MLLPT?branch=master)
<!-- badges: end -->

Maximum likelihood assignment of low-pass samples (MLLPT)

-----

#### Installation

``` r
# install.packages("devtools")
devtools::install_github('dombennett/treeman') # package requires newest version of treeman (not on CRAN)
devtools::install_github("T-Heide/MLLPT")
```

-----

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://t--heide.github.io/MLLPT/-informational)](https://t-heide.github.io/MLLPT/)

-----

#### Getting started

``` r
library(MLLPT)
library(phangorn)
```

##### 1) Input data

``` r
data("example_lp_data", package="MLLPT") # example dataset included in the package

# 1/3 a phylogenetic tree
print(example_lp_data$tree)
```

```text
Phylogenetic tree with 7 tips and 6 internal nodes.

Tip labels:
	EPICC_C518_A1_G4_D1, EPICC_C518_A1_G6_D1, EPICC_C518_B1_G6_D1, EPICC_C518_B1_G8_D1, EPICC_C518_D1_G6_D1, GL, ...

Rooted; includes branch lengths.
```


``` r
# 2/3 a phylo object (with id attribute)
print(example_lp_data$phydata)
```
```text
7 sequences with 7152 character and 19 different site patterns.
The states are 0 1 
```

``` r
head(attr(example_lp_data$phydata, "id")) # a non default attribute of mutation ids!
```
```text
[1] "1X0" "1X1" "1X2" "1X3" "1X4" "1X5"
```


``` r
# 3/3 a list of sample data
str(example_lp_data$samples)
```

```text
List of 2
  ..$ EPICC_C518_C1_G1_L1:'data.frame':	40452 obs. of  5 variables:
  .. ..$ alt_count : num [1:40452] 0 0 0 0 0 0 0 1 1 0 ...
  .. ..$ depth     : num [1:40452] 0 1 0 1 1 1 1 2 2 0 ...
  .. ..$ id        : chr [1:40452] "1X0" "1X1" "1X2" "1X3" ...
  .. ..$ cn_total  : num [1:40452] 2 2 2 2 2 2 2 2 2 2 ...
  .. ..$ cn_mutated: num [1:40452] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ EPICC_C518_D1_G1_L1:'data.frame':	36212 obs. of  4 variables:
  .. ..$ alt_count: num [1:36212] 0 1 0 1 0 0 0 1 0 1 ...
  .. ..$ depth    : num [1:36212] 0 1 1 1 3 0 0 2 1 4 ...
  .. ..$ id       : chr [1:36212] "1X0" "1X1" "1X2" "1X3" ...
  .. ..$ cn_total : num [1:36212] 2 2 2 2 2 2 2 2 2 2 ...
```

 
##### 2) Add LP samples

``` r
tree_with_lp_added = 
  with(example_lp_data,  {
    MLLPT::add_lowpass_sampled(
      tree = tree, 
      phydata = phydata, 
      sample_data = samples, 
      return_details = TRUE
    )
  })
```
 
```text
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=

Processing sample: EPICC_C518_C1_G1_L1 (1/2)


New values:

- Background rate: 0.01 -> 0 
- Purity: 1 -> 0.8313362 
- MLL: -349.1069 -> -192.6928 

=> Added sample (confidence: 1)

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=

Processing sample: EPICC_C518_D1_G1_L1 (2/2)


New values:

- Background rate: 0.01 -> 0.0007121943 
- Purity: 1 -> 0.6490007 
- MLL: -886.1944 -> -733.76 

=> Added sample (confidence: 1)

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=
```


##### 3) Plots

``` r
MLLPT::plot_tree(tree_with_lp_added)
```
<img src="man/figures/tree_plus_lp_example.png?raw=true" width="350">


``` r
MLLPT::plot_lp_loglik(tree_with_lp_added) + coord_flip()
```
<img src="man/figures/loglik_plus_lp_example.png?raw=true" width="200">


``` r
MLLPT::plot_lp_loglik_edge(tree_with_lp_added)
```
<img src="man/figures/loglik_edge_plus_lp_example.png?raw=true" width="850">


``` r
MLLPT::plot_sample_data(
  tree_with_lp_added, 
  external_purity_estimate = c(0.90, 0.48), 
  label_external_purity = "GT"
)
```
<img src="man/figures/sample_data_lp_example.png?raw=true" width="450">


-----

#### Copyright

Copyright (C) 2020, Timon Heide (timon.heide@icr.ac.uk)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.

-----

#### Contacts

Timon Heide, _Institute of Cancer Research, London, UK_.

[![](https://img.shields.io/badge/Email-timon.heide@icr.ac.uk-informational.svg?style=social)](mailto:timon.heide@icr.ac.uk)
[![](https://img.shields.io/badge/Github-T--Heide-informational.svg?style=social&logo=GitHub)](https://github.com/T-Heide)

