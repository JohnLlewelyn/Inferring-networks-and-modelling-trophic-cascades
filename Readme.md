## Naracoorte ecological network models
<img align="right" src="Naracoorte network.png" alt="Naracoorte Ecological Network" width="450" style="margin-top: 20px">

Scripts accompaning paper:

Llewelyn, J, G Strona, MC McDowell, CN Johnson, KJ Peters, DB Stouffer, SN de Visser, F Saltré, CJA Bradshaw. 2021. [Sahul’s megafauna species were vulnerable to plant-community changes due to their position in the trophic network](http://doi.org/10.1111/ecog.06089). <em>Ecography</em> doi:10.1111/ecog.06089

# Abstract
Extinctions stemming from environmental change often trigger trophic cascades and coextinctions. Bottom–up cascades occur when changes in the primary producers in a
network elicit flow-on effects to higher trophic levels. However, it remains unclear what determines a species’ vulnerability to bottom–up cascades and whether such cascades were a large contributor to the megafauna extinctions that swept across several continents in the Late Pleistocene. The pathways to megafauna extinctions are particularly unclear for Sahul (landmass comprising Australia and New Guinea), where extinctions happened earlier than on other continents. We investigated the potential role of bottom–up trophic cascades in the megafauna extinctions in Late Pleistocene Sahul by first developing synthetic networks that varied in topology to identify how network position (trophic level, diet breadth, basal connections) influences vulnerability to bottom–up cascades. We then constructed pre-extinction (~80 ka) network models of the ecological community of Naracoorte, south-eastern Sahul, to assess whether the observed megafauna extinctions could be explained by bottom–up cascades. Synthetic networks showed that node vulnerability to bottom–up cascades decreased with increasing trophic level, diet breadth and basal connections. Extinct species in the Naracoorte community were more vulnerable overall to these cascades than were species that survived. The position of extinct species in the network – tending to be of low trophic level and therefore having relatively narrow diet breadths and fewer connections to plants – made them vulnerable. However, these species also tended to have few or no predators, a network-position attribute that suggests they might have been particularly vulnerable to new predators. Together, these results suggest that trophic cascades and naivety to predators could have contributed to the megafauna extinction event in Sahul.

<strong>AUTHOR</strong>: John Llewelyn<br>
<strong>CONTACT</strong>: john.llewelyn@flinders.edu.au<br>
<strong>URL</strong>: <a href="http://GlobalEcologyFlinders.com">GlobalEcologyFlinders.com</a><br>
<strong>INSTITUTION</strong>: Flinders University<br>
<strong>RELEASE DATE</strong>: November 2021

## Guide to using scripts and data 

There are 2 main components:

1. Artificial networks. Includes models and script to calculate vulnerability to coextinction of nodes in these models (using Bayesian network and simulation approach). Use <code>script_toy_networks_vuln_edited.R</code> to run analyses calculating coextinction vulnerability.

2. Naracoorte network models. This includes scripts for:
- trophic niche space model development and validation (<code>Trophic niche space validation script.R</code>)
- building network models (<code>Build Naracoorte pre-extinction networks.R</code>)
- calculating coextincton vulnerability (<code>nara_vuln.py</code>)

Both components include a data folder containing the data needed to run the scripts.
