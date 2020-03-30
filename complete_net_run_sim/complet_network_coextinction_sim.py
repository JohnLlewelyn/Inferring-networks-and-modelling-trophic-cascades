import csv
from igraph import Graph
from numpy import inf
from random import sample,randrange,random
from collections import Counter
from scipy import trapz
from random import random,randrange,sample,shuffle
from numpy import arange,array,where,nan_to_num,zeros,mean,linspace
from copy import deepcopy
from scipy.spatial.distance import euclidean
from igraph import Graph
from numpy.random import normal,lognormal
from difflib import SequenceMatcher
import string
from numpy import percentile
from numpy import inf
import os as os

os.chdir('/Users/llew0024/Dropbox/trophic links/Docs/Final analyses/coextinction sim')

all_sp_18,all_sp_80 = set([]),set([])
for i in range(1,1001):
	f18 = [j for j in csv.reader(open('./data/nets_18/'+str(i)+'.csv','r'))] ###loop through file paths for each individual network
	f80 =[j for j in csv.reader(open('./data/nets_80_ok/'+str(i)+'.csv','r'))]
	for j in f18[1:]:
		all_sp_18 |= set(j[1:])
	for j in f80[1:]:
		all_sp_80 |= set(j[1:])
	print (i)


ext_spp = all_sp_80-all_sp_18   #list of extinct species

ext_dict=dict()                 #list of all species, indicating if they were extinct or not
for i in all_sp_80:
	if i in ext_spp:
		ext_dict[i] = 1
	else:
		ext_dict[i] = 0

####coextinction

def ext_casc(g,eee,pred_0_dict,prey_0_dict,bas,tre=0.70,rew=0.5):
	names = g.vs['name']                #g.vs from igraph - list of vertices (nodes) - attribute['names']
	try:
		gw_mat = g.get_adjacency(attribute='weight')	#igraph Graph object is transformed in a weighted adjacency matrix
	except:	#the try except here is lazy coding I used here to prevent the get_adjacency command to throw an error for an empty network as input (which happens in the last step of the co-extinction simulations); most likely not relevant for translating the code in R (but other similar issues might arise there).
		gw_mat = g.get_adjacency()
	gw_mat = array([[el for el in row] for row in gw_mat],dtype = float) #matrix is converted to array
	eee = [names.index(i) for i in eee]	#get the position in the array of the nodes going extinct
	new_mat=gw_mat.copy()	#create a copy of the original matrix
	pp = bas	#identify basal resources, as species having trophic level 0
	b_id=[node.index for node in g.vs if node['name'] in bas]
	pred_0 = array([pred_0_dict[i] for i in names]) #starting amount of resources using by consumers in the food web
	prey_0 = array([prey_0_dict[i] for i in names])	#starting amount of consumers using the resources in the food web
	ext=set(eee)	#transform array to set
	col_sum_0 = gw_mat.sum(0) #initial column totals; those correspond to the amount of resources used by consumers before the  co-extinction cascade starts
	while len(eee)>0:	#reiterate the steps in the co-extinction cascade as long as new secondary extinctions are triggered
		new_mat[eee,:]=0.0	#primary extinctions step 1: delete interactions where the extinct species act as consumers
		new_mat[:,eee]=0.0	#primary extinctions step 2: delete interactions where the extinct species act as resources
		col_sum_diff = new_mat.sum(0)/col_sum_0	#compute ratio between resources used by consumer after and before the co-extinction step
		col_sum_diff[where(col_sum_0==0)] = 1.0	#set to 1 the values in the ratio for which the initial value of used resources was 0
		new_mat = (col_sum_diff*new_mat.T).T	#the reduction in available resources for a given consumer is projected to the consumer's consumers, under the assumption that a reduction in available resources used by a consumer would result in a reduction of consumer's population
		#eee=set(where(new_mat.sum(0)<=tre)[0])-(ext|set(b_id))
		eee=set(where(new_mat.sum(0)/pred_0<=tre)[0])-(ext)	#identify secondary extinctions, as consumer losing a fraction of resources equal or larger than tre; this is evaluated on the basis of the resource used at the beginning of time, not at the beginning of the co-extinction cascade
		ext|=eee	#kept track of cumulated extinctions
		eee=array(list(eee))
		prey_avail=prey_0-new_mat.sum(1) #evaluate the amount of resources that have been freed by primary extinctions
		prey_avail[array(list(ext))]=0.0 #exclude from those the resources that have gone extinct
		pred_n = (new_mat>0).sum(1)	#compute the numbers of consumers using each resource in the network
		to_add=((prey_avail/pred_n)*rew)	#compute the amount of resources that can be reallocated to consumers; note that resources are reallocated only to species that are already using them, i.e. the structure of the network is not changed, but the weights of the consumers receiving the reallocable resource/s increase; allocable resource is equally divided between recipient consumers, and reduced by a 'rewiring' factor (rew).
		to_add[where(pred_n==0)]=0	#replace na values with 0
		to_add=(to_add*(new_mat.T>0))	#filter the to_add array to ensure the resource will be added to the proper set of consumers (i.e. those already using the resource, see previous comments)
		new_mat=(new_mat.T+to_add).T	#reallocate the resource
	g = Graph.Weighted_Adjacency([list(row) for row in gw_mat])	#convert the adjacency matrix back to an igraph Graph object
	g.vs['name'] = names
	g = g.subgraph(set(range(len(g.vs)))-ext)	#in the conversion process, nodes with 0 degree (i.e. with no in- or out-going links) are included in the network; this step takes them out
	ext_names = set([names[i] for i in ext]) #get the names (i.e. unique ids) of all extinct nodes/species
	b_id=[node.index for node in g.vs if node['name'] in bas] #identify basal resources (i.e. species in the food web having trophic level = 0)
	pre_check = set(g.vs['name'])	#get names of all nodes in the new network before the final step, where all nodes having no path to basal resources are deleted from the network
	to_keep=set([])
	for i in b_id:
		ddd=g.subcomponent(i, mode="OUT") #keep only species for which a path to basal resource exists
		to_keep|=set(ddd)
	g = g.subgraph(to_keep)	#reduce the graph to only species having a path connecting them to basal resources
	post_check = set(g.vs['name'])	#get the names of the reduced graph
	ext_names|=(pre_check-post_check) #add the species that have been removed by this last step to the set of extinct species
	return ext_names, g	#return the full set of species going extinct following the cascade, and the resulting food-web

#


###

out_seq = open('results/ext_seq_80_12_02_20.csv','a')
for net_n in range(1,11):
	tre = random()
	f80 = [j for j in csv.reader(open('./data/nets_80_ok/'+str(net_n)+'.csv','r'))]
	net = Graph.TupleList([i[::-1] for i in f80[1:]],directed=True)
	p_id = net.vs['name'].index('plants')
	i_id = net.vs['name'].index('invertebrates')
	f_id = net.vs['name'].index('FishInv')
	p_n,i_n,f_n = [],[],[]
	net = []
	vert = []
	for i in f80[1:]:
		res,cons = i[1:][::-1]
		net.append([res,cons])
		if res == 'plants':
			p_n.append(cons)
		if res == 'invertebrates':
			i_n.append(cons)
		if res == 'FishInv':
			f_n.append(cons)
		if res not in ['plants','invertebrates','FishInv']:
			vert.append(res)
		if cons not in ['plants','invertebrates','FishInv']:
			vert.append(cons)
	vert = set(vert)
	co_p = Counter(p_n)
	co_i = Counter(i_n)
	co_f = Counter(f_n)
	min_p,max_p = max(co_p.values()),sum(co_p.values())
	min_i,max_i = max(co_i.values()),sum(co_i.values())
	min_f,max_f = max(co_f.values()),sum(co_f.values())
	p_set = ['p'+str(i) for i in range(int(min_p*(1+random())))]#int(max_p*ov))]
	i_set = ['i'+str(i) for i in range(int(min_i*(1+random())))]#int(max_i*ov))]
	f_set = ['f'+str(i) for i in range(int(min_f*(1+random())))]#int(max_f*ov))]
	comp_net = []
	for i in net:
		res,cons = i
		if res == 'plants':
			res = sample(p_set,1)[0]
		if res == 'invertebrates':# and random()<0.75:
			res = sample(i_set,1)[0]
		if res == 'FishInv':
			res = sample(f_set,1)[0]
		comp_net.append([res,cons])
	non_vert = set(p_set+f_set+i_set)
	conn_ip = random()/10.0
	for i in range(int(len(p_set)*len(i_set)*conn_ip)):
		comp_net.append([sample(p_set,1)[0],sample(i_set,1)[0]])
	for i in range(int(len(p_set)*len(f_set)*conn_ip)):
		comp_net.append([sample(p_set,1)[0],sample(f_set,1)[0]])
	for i in i_set:
		if i not in [j[1] for j in comp_net]:
			comp_net.append([sample(p_set,1)[0],i])
	for i in f_set:
		if i not in [j[1] for j in comp_net]:
			comp_net.append([sample(p_set,1)[0],i])
	lll = array([lognormal() for i in range(500000)]) #generate 10000 numbers sampled from a lognormal distribution
	lll/=max(lll) #get the max value of the lognormal
	g = Graph.TupleList(comp_net,directed=True)
	g = g.simplify()	#eliminate duplicate edges
	g.es['weight'] = sample(list(lll),len(g.es))
	exp_name = ['random']#,'best','worst'] #names of the three scenarios
	for exp in range(1): #iterate over the scenarios
		net = g.copy() #create a copy of the network under study
		wg = net.get_adjacency(attribute='weight')	#convert food web to weighted adjacency matrix
		wg = array([[cell for cell in row] for row in wg],dtype = float)	#convert matrix to array
		names = net.vs['name']	#get names from network
		spn = len(names) #number of nodes/species
		pred_sum = wg.sum(0) #get sums of columns, i.e. total amounts of resources used by consumers
		prey_sum = wg.sum(1)	#get sums of rows, i.e. total amounts of consumers using a resource
		pred_0_dict = dict([[names[i],pred_sum[i]] for i in range(spn)])	#make dictionary mapping values to species ids
		prey_0_dict = dict([[names[i],prey_sum[i]] for i in range(spn)])	#make dictionary mapping values to species ids
		bas = set(array(net.vs['name'])[where(array(net.degree(mode='IN'))==0)]) #identify basal resources, as nodes with only outgoing links but no incoming ones;
		pe = 0.0 #primary extinction counter
		for ext_ev in range(len(g.vs)): #iterate over the entire length of the network
			if len(net.vs)>0:	#evaluate co-extinction only if some node is still in the network
				eee = sample(net.vs['name'],1) #in the random scenario, just sample a random node
				ext_names, net = ext_casc(net,eee,pred_0_dict,prey_0_dict,bas,tre=tre,rew=0.5) #evaluate co-extinctions stemming from the primary extinction (the single node in the eee list); the function return the complete list of nodes that went extinct (including the node in eee); and the network after the co-extinction cascade (ready for the next node-removal/co-extinction step).
				co_exts = ext_names-(set(eee)|non_vert) #get the names of nodes that went extinct following the co-extinction (i.e. take all the extinct nodes indicated by the function ext_casc, and remove the node that we removed, i.e. the one in list eee
				if co_exts!=set([]): #if co-extinction happens, write those, one per line, to the output file, indicating the network id (rep), the scenario (exp_name[exp]), the co-extinct species (co_ext_sp), and after how many node removal step the co-extinction happened (ext_ev)
					for co_ext_sp in co_exts:
						out_seq.write(','.join(map(str,[net_n,ext_dict[co_ext_sp],co_ext_sp,ext_ev]))+'\n')
			pe+=1 #update primary extinction counter
	print (net_n) #print network number to track progress in the analysis


out_seq.close()


######################NET PROPS
out = open('for_mds.csv','w')
out.close()
for net_n in range(1,11):
	f80 = [j for j in csv.reader(open('./nets_80_ok/'+str(net_n)+'.csv','r'))]
	net = Graph.TupleList([i[::-1] for i in f80[1:]],directed=True)
	p_id = net.vs['name'].index('plants')
	i_id = net.vs['name'].index('invertebrates')
	f_id = net.vs['name'].index('FishInv')
	p_n,i_n,f_n = [],[],[]
	net = []
	vert = []
	for i in f80[1:]:
		res,cons = i[1:][::-1]
		net.append([res,cons])
		if res == 'plants':
			p_n.append(cons)
		if res == 'invertebrates':
			i_n.append(cons)
		if res == 'FishInv':
			f_n.append(cons)
		if res not in ['plants','invertebrates','FishInv']:
			vert.append(res)
		if cons not in ['plants','invertebrates','FishInv']:
			vert.append(cons)
	vert = set(vert)
	co_p = Counter(p_n)
	co_i = Counter(i_n)
	co_f = Counter(f_n)
	min_p,max_p = max(co_p.values()),sum(co_p.values())
	min_i,max_i = max(co_i.values()),sum(co_i.values())
	min_f,max_f = max(co_f.values()),sum(co_f.values())
	p_set = ['p'+str(i) for i in range(int(min_p*(1+random())))]#int(max_p*ov))]
	i_set = ['i'+str(i) for i in range(int(min_i*(1+random())))]#int(max_i*ov))]
	f_set = ['f'+str(i) for i in range(int(min_f*(1+random())))]#int(max_f*ov))]
	comp_net = []
	for i in net:
		res,cons = i
		if res == 'plants':
			res = sample(p_set,1)[0]
		if res == 'invertebrates':# and random()<0.75:
			res = sample(i_set,1)[0]
		if res == 'FishInv':
			res = sample(f_set,1)[0]
		comp_net.append([res,cons])
	non_vert = set(p_set+f_set+i_set)
	conn_ip = random()/10.0
	for i in range(int(len(p_set)*len(i_set)*conn_ip)):
		comp_net.append([sample(p_set,1)[0],sample(i_set,1)[0]])
	for i in range(int(len(p_set)*len(f_set)*conn_ip)):
		comp_net.append([sample(p_set,1)[0],sample(f_set,1)[0]])
	for i in i_set:
		if i not in [j[1] for j in comp_net]:
			comp_net.append([sample(p_set,1)[0],i])
	for i in f_set:
		if i not in [j[1] for j in comp_net]:
			comp_net.append([sample(p_set,1)[0],i])
	net = Graph.TupleList(comp_net,directed=True)
	net = net.simplify()
	g = net.copy()
	ext_codes = [i for i in range(len(g.vs)) if g.vs['name'][i] in ext_spp]
	surv_codes = [i for i in range(len(g.vs)) if g.vs['name'][i] in vert-ext_spp]
	pr = g.pagerank(directed=True)
	bet = g.betweenness(directed=True)
	eig = g.eigenvector_centrality(directed=True)
	close_in = g.closeness(mode='IN')
	core_in = g.coreness(mode='IN')
	deg_in = g.degree(mode='IN')
	ecc_in = g.eccentricity(mode='IN')
	close_out = g.closeness(mode='OUT')
	core_out = g.coreness(mode='OUT')
	deg_out = g.degree(mode='OUT')
	ecc_out = g.eccentricity(mode='OUT')
	out = open('for_mds.csv','a')
	for i in ext_codes+surv_codes:
		status = 'surv'
		if i in ext_codes:
			status = 'ext'
		out.write(','.join(map(str,[
				status,
				g.vs['name'][i],
				pr[i],
				bet[i],
				eig[i],
				close_in[i],
				core_in[i],
				deg_in[i],
				ecc_in[i],
				close_out[i],
				core_out[i],
				deg_out[i],
				ecc_out[i]]))+'\n')
	out.close()
	print (net_n)
