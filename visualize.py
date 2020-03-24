import importlib, sys

fname = sys.argv[1]
sys.argv[1] = 'Weibull'	# to fix evalulate.py error

from evaluate import *
from matplotlib import pyplot as plt
from numpy import inf

colors = {
	'bat':		'#020202',
	'bee':		'#bbaf22',
	'cuckoo':	'#ec994a',
	'firefly': 	'#ff0000',
	'fish':		'#3d7bff',
	'pollination':'#ff3dff',
	'pso':		'#3dd6ff',
	'wolf':		'#e03e3e',
	'NONE':		'#00aa00'
}



def decode(bitstring):
	'''
	Takes some bit-pattern (some member of the population),
	gives back the corresponding algorithms (phase 1, 2, 3)
	as well as the parameters passed to ph1
	'''
	groups = [bitstring[x*nsga_bpp:][:nsga_bpp] for x in range( int( len(bitstring) / nsga_bpp ) )]

	index_1 = round(bin_signum(groups[0], (len(stage1)-1)/2, (len(stage1)-1)/2))
	index_2 = round(bin_signum(groups[1], (len(stage2)-1)/2, (len(stage2)-1)/2))
	index_3 = round(bin_signum(groups[2], (len(stage3)-1)/2, (len(stage3)-1)/2))
											# get indexes of each stage, by rounding
											# stage_num / 2 +- stage_num / 2

	pop_size = round(bin_signum(groups[3], model_pop_count[0], model_pop_count[1]))
	gen_count = round(bin_signum(groups[4], model_generations[0], model_generations[1]))
											# get bits corresponding to model gens and size
											# and generate pop size / generations based on values

	#model_search_space = [model["search_space"] for i in range(model["dimensions"])]

	params = (lambda x:x , (), gen_count, [0]*pop_size)
											# get parameters for swarm algorithms
	param_bits = groups[5:]
	for idx, prm in enumerate(stage2[index_2]["params"]):
		this_bits = param_bits[idx]
											# get part of bitstring corresponding to parameter
		if nsga_free_range:
			params += (bin_signum(this_bits, 0.5, 0.5),)
		else:
			params += (bin_signum(this_bits, prm[0], prm[1]),)

	return stage1[index_1], stage2[index_2]["algo"], stage3[index_3], params	# send back algorithms & parameters for timing



pops = []
			# take pops from file and import to script
with open(fname) as f:
	for gen_line in f:
		pops.append([])
		linedata = gen_line.split(" ")
		candidate_count = len(linedata) / 3	#groups of 3
		for i in range(int(candidate_count)):
			pops[-1].append( {'bitstring': linedata[i*3] , 'objectives': [float(linedata[i*3 + 1]), float(linedata[i*3 + 2])] } )

allpops = []
for g in pops:
	for p in g:
		allpops.append(p)



pop = pops[-1]

pop.sort(key = lambda x: (x["objectives"][0], x["objectives"][1]))	# sort by conv then time

sep = "\t"

print('\033[4m' + sep.join(['runtime', 'avg conv', 'algo', 'conv method', 'gens', 'pop', 'algo params', 'init estimate score']) + "\\\\" + '\033[0m')
for p in allpops:
	r, r2, r3, params = decode(p["bitstring"])
	n = r2.__module__ if r2 != None else "NONE"
	c = colors[n]

	#print(sep.join([str(round(x,12)).ljust(12, ' ') for x in p["objectives"]]), end=sep)

	#print(f"{n[:6] if r2 != None else 'NONE'}{sep}{r3.__name__[:6] }{sep}{params[2] if r2 != None else 'n/a'}{sep}{len(params[3])if r2 != None else 'n/a'}", end=sep)
	if r2 != None:
		for i in params[4:]:
			pass
			#print(round(i,3), end=sep)

	#print(p["bitstring"], end=sep)
	print(p['objectives'][0], p['objectives'][1])
	#print("\\\\" if sep != "\t" else "")





'''				# plot pop pcts
algs = []
for gen, pop in enumerate(pops):
	t = {'NONE':0}								# code to visualize population percentages of each algorithm
	for i in [ x[:-3] for x in os.listdir('algorithms') if x[-3:] == '.py']:
		t[i] = 0
	for p in pop:
		_,b,_,_ = decode(p["bitstring"])
		m = b.__module__ if b != None else 'NONE'
		t[m] += 1
	algs.append(t)

names = [x[:-3] for x in os.listdir('algorithms') if x[-3:] == '.py']
names.append('NONE')
used = []
prog = {}
for n in names:
	prog[n] = [x[n] for x in algs]
for n in names:
	t = n
	if not(t == 'bee' or t == 'NONE'):
		t = 'Other'
	plt.step([l/nsga_pop_size for l in prog[n]],'#000000',label=t.lower() if t not in used else None, ls='-' if n=='NONE' else ':' if n=='bee' else '--')
	used.append(t)


plt.title('Algorithm Population Share Ratio')
plt.legend(loc='best',
#plt.xticks([0] + [x-1 for x in range(16,129,16)], ['1'] + [str(x) for x in range(16,129,16)])
plt.ylabel('Percentage of Population')
plt.xlabel('Generation')
plt.show()
'''


plt.ticklabel_format(useOffset=False)
buckets = []
dsets = ['SYS1','SYS2','SYS3','S2',	'S27','SS3','CSR1','CSR2','CSR3']
marks = ['.',	'^',	'1',  's',	'P',  '*',  'x',   'd',   '>']

colors = ['0','0.65','0.65']

for i in range(9):
	p0 = allpops[i*3 +0]['objectives']
	p1 = allpops[i*3 +1]['objectives']
	p2 = allpops[i*3 +2]['objectives']
	buckets.append(( (p0[0] + p1[0] + p2[0])/3, (p0[1] + p1[1] + p2[1])/3 ))

if False:	# plot avg or all

	for i, p in enumerate(buckets):
		if i == 0:
			continue
		plt.plot(p[0], p[1], marker=marks[i], color='k', label=dsets[i])

else:
	#for p in pops:
	for ind, p in enumerate(allpops):
		if ind == 0:
			continue
		set_ind = int(ind/3)
		dset = dsets[set_ind]
		color = colors[ind%3]
		ecolor = None #if color == 'k' else 'k'
		plt.scatter(p['objectives'][0], p['objectives'][1], marker=marks[set_ind], color=color, edgecolors=ecolor, label=dsets[set_ind] if ind%3==0 else None)
		plt.ylim([1-0.000125, 1.002125])
		
#plt.xlim([0,0.015])
plt.xlabel('Runtime (s)',fontsize=12)
plt.ylabel('Accuracy ($1+\epsilon$)',fontsize=12)

plt.legend(loc='best',scatterpoints = 1)
plt.show()





