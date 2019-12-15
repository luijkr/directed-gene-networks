# Repository for manuscript _"Genome-wide identification of directed gene networks using large-scale population genomics data"_ by Luijk _et al_.

#### Construct genetic instruments

_construct_genetic_instruments.R_ loads the data and creates genetic instruments for all genes with nearby genetic variants in our data. The lasso implemented in _glmnet_ is used.

#### Identify _trans_-effects

_trans_effects.R_ identifies _trans_-effects. This script is based on work bij Sikorska _et al_. (https://doi.org/10.1186/1471-2105-14-166).

#### Simulations

_simulation.R_ runs the simulation study. Each combination of settings is run 500 times.