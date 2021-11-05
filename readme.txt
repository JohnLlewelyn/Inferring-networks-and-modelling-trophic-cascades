## Naracoorte ecological network models

Scripts accompaning paper:

Llewelyn, J, G Strona, MC McDowell, CN Johnson, KJ Peters, DB Stouffer, SN de Visser, F Saltré, CJA Bradshaw. 2021. <a href="10.1111/ecog.06089">Sahul’s megafauna species were vulnerable to plant-community changes due to their position in the trophic network</a>. <em>Ecography</em> doi:10.1111/ecog.06089

## Guide to using script and data 

There are 2 main components:

1. Artificial networks. Includes models and script to calculate vulnerability to coextinction of nodes in these models (using Bayesian network and simulation approach). Use script_toy_networks_vuln_edited.R to run analyses calculating coextinction vulnerability.

2. Naracoorte network models. This includes scripts for:
- trophic niche space model development and validation (Trophic niche space validation script.R)
- building network models (Build Naracoorte pre-extinction networks.R)
- calculating coextincton vulnerability (nara_vuln.py)

Both components include a data folder containing the data needed to run the scripts.
