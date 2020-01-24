import importlib, sys
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

def decode(bitstring, model = models.models[sys.argv[1]]):
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

	model_search_space = [model["search_space"] for i in range(model["dimensions"])]

	params = (model["objective"], model_search_space, gen_count, [0]*pop_size)
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
with open(sys.argv[2]) as f:
	for gen_line in f:
		pops.append([])
		linedata = gen_line.split(" ")
		candidate_count = len(linedata) / 3	#groups of 3
		for i in range(int(candidate_count)):
			pops[-1].append( {'bitstring': linedata[i*3] , 'objectives': [float(linedata[i*3 + 1]), float(linedata[i*3 + 2])] } )


pop = pops[-1]

pop.sort(key = lambda x: (x["objectives"][0], x["objectives"][1]))	# sort by conv then time

sep = "\t"

print('\033[4m' + sep.join(['runtime', 'avg conv', 'algo', 'conv method', 'gens', 'pop', 'algo params', 'init estimate score']) + "\\\\" + '\033[0m')
for p in pop:
	r, r2, r3, params = decode(p["bitstring"])
	n = r2.__module__ if r2 != None else "NONE"
	c = colors[n]

	print(sep.join([str(round(x,12)).ljust(12, ' ') for x in p["objectives"]]), end=sep)
	print(f"{n[:6] if r2 != None else 'NONE'}{sep}{r3.__name__[:6] }{sep}{params[2] if r2 != None else 'n/a'}{sep}{len(params[3])if r2 != None else 'n/a'}", end=sep)

	if r2 != None:
		for i in params[4:]:
			print(round(i,3), end=sep)
	print("\\\\" if sep != "\t" else "")




'''
types = {}		# bar chart of none types
for p in pop:
	ptype, conv = decode(p['bitstring'])[1:][:2]
	if ptype == None:
		conv = conv.__name__
		if conv in types:
			types[conv] += 1
		else:
			types[conv] = 1
plt.bar(list(range(len(types))),
		[types[x] for x in types],
		tick_label = [str(x) for x in types.keys()])	
plt.title('Subset of NONE types')
plt.show()
'''

'''
for idx, pop in enumerate(pops):	# plot runtime changes per generation
	for p in pop:
		plt.plot(idx+2, p['objectives'][0], 'r*')
	plt.plot(idx+2, sum([x['objectives'][0] for x in pop])/nsga_pop_size, 'b')
plt.title('All runtimes versus generation')
plt.show()
'''

'''#								# plot best candidates of all runs
names = [x[:-3] for x in os.listdir('algorithms') if x[-3:] == '.py']
names.append('NONE')
prog = {}
for n in names:
	prog[n] = 0
for p in frontrunners:
	_, s2, _, _ = decode(p["bitstring"])
	s2 = s2.__module__ if s2 != None else "NONE"
	if s2 != None:
		prog[s2] += 1
	else:
		prog['NONE'] += 1
klist = list(prog.keys()).copy()
for t in klist:
	if prog[t] == 0:
		del prog[t]
keys = prog.keys()
values = [prog[x] for x in keys]
plt.bar(list(range(len(keys))), values, align='center', alpha = 0.5)
plt.xticks(list(range(len(keys))), keys)
plt.ylabel('Frontrunner algorithm type')
plt.show()
'''


'''# 			PLOT SIZE AND GENERATIONS VS ITERATIONS
plt.plot([sum([x[0] for x in l])/len(l) for l in pops] , 'b', label='pop size')
plt.plot([sum([x[0] for x in l])/len(l) for l in gens], 'r', label='gen count')
plt.title('Population Size & Swarm Generations vs Iterations (Average)')
plt.legend(loc='best')
plt.show()
'''

'''
				#PLOT TIME/EPSILON VS ITERATIONS
ers = []
tss = []
for iteration, pt in enumerate(pops):
	er = 0
	ts = 0
	c = 1
	for memb in pt:
		plt.plot(iteration, memb['objectives'][0],'r*')
		ts += memb['objectives'][0]
		if (not np.isnan(memb['objectives'][1])) and (0 < memb['objectives'][1] < 3):
			plt.plot(iteration, memb['objectives'][1]-1,'b*')
			er += memb['objectives'][1]-1
			c += 1
	er /= c
	ts /= len(pt)
	ers.append(er)
	tss.append(ts)
plt.plot(tss, 'r',label='time')
plt.plot(ers, 'b', label='epsilon')
plt.title('Time & Error vs Iterations')
plt.legend(loc='best')
plt.show()
'''

				# plot pop pcts
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
plt.legend(loc='best')
plt.xticks([0] + [x-1 for x in range(16,129,16)], ['1'] + [str(x) for x in range(16,129,16)])
plt.ylabel('Algorithm Share Percentage')
plt.xlabel('Generation')
plt.show()



names = []
for p in pops[-1]:
	r1, r2, r3, prop = decode(p['bitstring'])
	n = r2.__module__ if r2 != None else "NONE"


	c = 'k'
	if n == 'NONE':
		c = '.'
	if n == 'bee':
		c = 'x'
	if str(round(p['objectives'][1],12)) == '1.000325607232':
		c = '+'
	plt.plot(p['objectives'][0], p['objectives'][1], marker=c, color='k', label=None if n in names else n.lower())
	if n not in names:
		names.append(n)

#plt.ylim([0.999, 1.01])
#plt.xlim([0,0.025])
plt.title('Error & Runtime Pareto Plot')
plt.xlabel('Runtime (s)')
plt.ylabel('1+epsilon error')
plt.legend(loc='best')
plt.show()
	
'''	ANIMATE PARETO PLOT, DONT DO THIS
import matplotlib
matplotlib.use('Agg')
import imageio
ims = []
for i, pop in enumerate(pops):

	fig, ax = plt.subplots()
					# plot pareto curve
	for p in pop:
		ax.plot(p['objectives'][0], p['objectives'][1], 'b*')

	#plt.ylim([0.99, 1.1])
	#plt.xlim([0,0.025])
	#ax.title('Error vs Runtime Pareto Front')
	#ax.xlabel('Runtime (s)')
	#ax.ylabel('1+epsilon error')
	ax.set(xlabel='Runtime (s)', ylabel='1+epsilon error', title=f'Error vs Runtime Pareto Curve: {i+1}')
	
	ax.set_ylim(0, 8)
	ax.set_xlim(0,0.04)
	fig.canvas.draw()
	img = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
	img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))

	ims.append(img)
	print(f'{i}  ',end="\r")

kwargs_write = {'fps':16, 'quantizer':'nq'}
imageio.mimsave('./powers.gif', ims, fps=16)
'''


'''
front0 = []
front1 = []
for iteration, popu in enumerate(pops):
	front = min(popu, key=lambda x: x['objectives'][0]*x['objectives'][1])
	front0.append(front['objectives'][0]*30)
	front1.append(front['objectives'][1]-1)

plt.plot(front0, 'b', label='runtime')
plt.plot(front1, 'r', label='error')

plt.title('Frontrunner accuracy and runtime')
plt.legend()
plt.show()
'''