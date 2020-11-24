Guide to using script and data 
	
There are five main components:

1) trophic niche space model development and validation
2) inferring network models of the Late Pleistocene Naracoorte assemblage
3) adding information on plant, invertebrate, and aquatic animal diversity, and then adding links to and from these nodes
4) evaluating vulnerability to bottom-up coextinction cascades using an algorithm that iteratively calculates the product of the vulnerabilities of the focal  	 node’s resource nodes
5) assess the network position of extinct versus surviving vertebrate nodes

1. Trophic niche space model development and validation
Uses GloBI and Serengeti data to define trophic niche space models and measure how much these models over-estimate number of prey so adjustments can be made.

Script (R): 
- Trophic niche space validation script.R

Data
- interactions_globi_cat,csv
- data in training folder (from GloBI)
- data in validation folder (from GloBI)
- Predator-Prey Data.csv
- ll.prey.csv #or can generate this data using rglobi
- Baskerville et al Serengeti.csv
- Baskerville et al species.csv
- Baskerville animals.csv
- MamFuncDat.csv
- GLOBI.BS.prey.csv

2. Inferring network models of the Late Pleistocene Naracoorte assemblage
Uses the best trophic niche space model identified in 1) and the Late Pleistocene Naracoorte assemblage to build network models.

Script (R):
- Build Naracoorte pre-extinction networks.R
- megafauna vetting.R (sourced in previous script)

Data:
- keep80.rds
- Serengeti.density.rds (generated in Trophic niche space validation script.R)
- prey.ob.rds (generated in Trophic niche space validation script.R)
- mammal.herb.density.rds (generated in Trophic niche space validation script.R)
- Table S5 diet breadth.csv

3. and 4. Convert RDS network files to csv files using R. Then with python, add information on plant, invertebrate and aquatic animal diversity, and use an algorithm to estimate vulnerability to bottom-up cascades.
Script (R):
- extract.R

Script (Python):
- co_extinctions_20_11_20.py

Data:
- x80nets.rds (generated in previous step) 
- x18lump.rds (post-extinction networks, attached)

5. Compare position of extinct versus surviving species in the network

Script (R):
- pca_position.R

Data:
- network_position_metrics.rds(from complete networks)
