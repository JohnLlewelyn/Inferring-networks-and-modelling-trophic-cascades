## Naracoorte ecological network models

Scripts accompaning paper:

Llewelyn, J, G Strona, MC McDowell, CN Johnson, KJ Peters, DB Stouffer, SN de Visser, F Saltré, CJA Bradshaw. 2021. [Sahul’s megafauna species were vulnerable to plant-community changes due to their position in the trophic network](http://doi.org/10.1111/ecog.06089). <em>Ecography</em> doi:10.1111/ecog.06089

## Guide to using scripts and data 

There are 2 main components:

1. Artificial networks. Includes models and script to calculate vulnerability to coextinction of nodes in these models (using Bayesian network and simulation approach). Use <code>script_toy_networks_vuln_edited.R</code> to run analyses calculating coextinction vulnerability.

2. Naracoorte network models. This includes scripts for:
- trophic niche space model development and validation (<code>Trophic niche space validation script.R</code>)
- building network models (<code>Build Naracoorte pre-extinction networks.R</code>)
- calculating coextincton vulnerability (<code>nara_vuln.py</code>)

Both components include a data folder containing the data needed to run the scripts.
