#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 09:29:12 2020

@author: giovanni strona
"""
import csv
from igraph import Graph
from random import sample,random
from numpy import array,where,inf,mean
from numpy.random import pareto
from copy import deepcopy

all_sp_18,all_sp_80 = set([]),set([])
for i in range(1,1001):
	f18 = [j for j in csv.reader(open('./nets_18/'+str(i)+'.csv','r'))]
	f80 =[j for j in csv.reader(open('./nets_80/'+str(i)+'.csv','r'))]
	for j in f18[1:]:
		all_sp_18 |= set(j[1:])
	for j in f80[1:]:
		all_sp_80 |= set(j[1:])
	print (i)


ext_spp = all_sp_80-all_sp_18

ext_dict=dict()
for i in all_sp_80:
	if i in ext_spp:
		ext_dict[i] = 1
	else:
		ext_dict[i] = 0



net_n = 1
f80 = [j for j in csv.reader(open('./nets_80/'+str(net_n)+'.csv','r'))]
net = [i[1:][::-1] for i in f80[1:]]
g = Graph.TupleList(net,directed=True)
spp = g.vs['name']



all_paths = g.get_shortest_paths('plants',spp,mode='OUT')

tl = []
for t_sp in range(len(g.vs)):
	#tp_all = g.shortest_paths_dijkstra(p_set_ok,t_sp)
	row = []
	for p_sp in ['plants']:
		rg = g.copy()
		tp = rg.get_shortest_paths(p_sp,t_sp)[0]
		while len(tp)>1:
			ttl = len(tp)
			row.append(ttl)
			rg.delete_edges([tuple(tp[-2:])])
			tp = rg.get_shortest_paths(p_sp,t_sp,mode='OUT')[0]
	if row!=[]:
		tl.append([mean(row),min(row),max(row)])
	else:
		tl.append([1,1,1])
	print (g.vs['name'][t_sp],tl[-1])



tl_dict = dict([[g.vs['name'][i],tl[i]] for i in range(len(g.vs))])


par_dist = [int(round(i)+1) for i in (pareto(a=1.02,size=100000)) if i<=51]

####coextinction
out = open('results_22_07_21.csv','w')
out.write('net_n,extinct_status,species,vulnerability,plant_n,mean_tl,min_tl,max_tl,paths_to_bas,deg_in,deg_out\n')
for net_n in range(1,1001):
	f80 = [j for j in csv.reader(open('./nets_80/'+str(net_n)+'.csv','r'))]
	net = [i[1:][::-1] for i in f80[1:]]
	p_set = ['p'+str(i) for i in range(1300)]
	i_set = ['i'+str(i) for i in range(6000)]
	f_set = ['f'+str(i) for i in range(23)]
	#create links between vertebrates and plants/invertebrates
	comp_net = []
	for i in net:
		res,cons = i
		if cons not in ['plants','invertebrates','fish']:
			if res == 'plants':
				res = sample(p_set,1)[0]
			if res == 'invertebrates':# and random()<0.75:
				res = sample(i_set,1)[0]
			if res == 'fish':
				res = sample(f_set,1)[0]
		comp_net.append([res,cons])
	#create links between plants and invertebrates
	non_vert = set(p_set+f_set+i_set)
	for i_sp in i_set:
		for j in range(sample(par_dist,1)[0]):
			comp_net.append([sample(p_set,1)[0],i_sp])
	for f_sp in f_set:
		for j in range(sample(par_dist,1)[0]):
			comp_net.append([sample(p_set,1)[0],f_sp])
	pres_nodes = set([j[1] for j in comp_net])
	for i in set(i_set):
		if i not in pres_nodes:
			comp_net.append([sample(p_set,1)[0],i])
	for i in set(f_set):
		if i not in [j[1] for j in comp_net]:
			comp_net.append([sample(p_set,1)[0],i])
	all_non_an = set(p_set)|set(i_set)|set(f_set)
	g = Graph.TupleList(comp_net,directed=True)
	g.delete_vertices(['plants','invertebrates','fish'])
	plants = sorted(list(set(p_set)&set(g.vs['name'])))
	non_plants = [i for i in g.vs['name'] if i not in plants]
	ddd = g.shortest_paths_dijkstra(non_plants,plants, mode="IN")
	to_keep = [non_plants[i] for i in range(len(non_plants)) if min(ddd[i])!=inf]+plants
	g = g.subgraph(to_keep)	#reduce the graph to only species having a path connecting them to basal resources
	g = g.simplify()	#eliminate duplicate edges
	ddd = g.shortest_paths_dijkstra(target=plants, mode="IN")
	deg_in = g.degree(mode='IN')
	deg_out = g.degree(mode='OUT')
	deg_dict = dict([[g.vs['name'][k],[len([j for j in ddd[k] if j!=inf]),deg_in[k],deg_out[k]]] for k in range(len(g.vs))])
	spp = g.vs['name']
	out_net = open('./networks/'+str(net_n)+'.csv','w')
	for i in g.es:
		v1,v2 = i.tuple
		out_net.write(spp[v1]+','+spp[v2]+'\n')
	out_net.close()
	plants = sorted(list(set(p_set)&set(g.vs['name'])))
	an_spp = set(spp)-all_non_an
	g_ = g.copy()
	sc=1
	plants_ = deepcopy(plants)
	plant_n = len(plants)
	for p in sample(plants,plant_n):
		spp_ = g_.vs['name']
		an_spp_ = set(spp_)-all_non_an
		g_.delete_vertices(p)
		plants_.remove(p)
		if len(plants_)>0:
			ddd = g_.shortest_paths_dijkstra(target=plants_, mode="IN")
			to_del = [j for j in range(len(ddd)) if min(ddd[j])==inf]
		else:
			to_del = g_.vs['name']
		g_.delete_vertices(to_del)
		coe = an_spp_-set(g_.vs['name'])
		for coe_sp in coe:
			out.write(','.join(map(str,[net_n,ext_dict[coe_sp],coe_sp,sc,plant_n]+tl_dict[coe_sp]+deg_dict[coe_sp]))+'\n')
		sc+=1
	print (net_n) #print network number to track progress in the analysis


out.close()




