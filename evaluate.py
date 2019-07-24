# TODO - histogram of all fn runtimes, "hills" should approach left at some point
# TODO - check if converged uses endpoints of ranges (perhaps to increase ranges)
# TODO - explore variants of best ones (bee, bat, pso, etc)
# TODO - separate each model in file into model folder

# TODO - add bits for skipping swarm algo
# TODO - use both CV datasets
# TODO - add estimate fn for wei (a = n, b = n/sum(bi), c = 1)
# TODO - check probability of range selection, make sure it's random
# TODO - write nsga pgh of running with params
# TODO - fix table in paper (one name per algo, zfill ranges, re-order to match paper)
# TODO - before illustrations, talk about methods then covariates
# TODO - table of NSGA variables
# TODO - paragraph on open nsga constraints

# TODO - 1+ep and time vs iterations graph
# TODO - pop and gens vs iterations (swarm)
# TODO - add no swarm to percentage graph

import matplotlib.pyplot as plt, csv
import sys, os, models, nsga, numpy as np

sys.path.insert(0, f"algorithms")
import bat, bee, cuckoo, firefly, fish, pollination, pso, wolf

#---- NSGA SETTINGS ------------------------------

nsga_max_gens = 32				# number of generations in NSGA-II
nsga_pop_size = 64 				# how many combinations there are, must be even
nsga_p_cross = 0.98				# mutation crossover probability
nsga_fn_evals = 16				# how many evaluations to average
nsga_bpp = 16 					# bits / precision to use for each parameter
nsga_cross_breed = False		# allow algorithms to change during the process
nsga_free_range = True			# sets all parameter ranges to 0.5 +- 0.5 instead of predefined (needs more gens)
nsga_cross_head_tail = True		# uses head/tail crossover instead of random bits

model_pop_count = [16, 10]		# pop size for method
model_generations = [12, 10]	# generations used in method


stage1 = [lambda x: x]*2
stage2 = [	{	
				"algo":	firefly.search,
				"params":[
					[0.95, 0.05],	[0.96, 0.04]
				]
			},
			{						# algo points to search function in swarm file
				"algo":	pso.search,	# contains all constant parameters in terms of 
				"params":[			# [A, B] where the passed value is A + B*f(bits)
					[0.5, 0.15],	[0.1, 0.05],	[0.1, 0.05]
				]
			},
			{
				"algo": cuckoo.search,
				"params":[
					[0.97, 0.03],	[0.97, 0.03],	[0.25, 0.1]
				]
			},
			{
				"algo": bat.search,
				"params":[
					[0.1, 0.1],		[0.9, 0.1],		[0.9, 0.08],	[0.8, 0.15]
				]
			},
			{
				"algo": bee.search,
				"params":[
					[0.33, 0.10],	[0.33, 0.10],	[0.50, 0.25],	[0.50, 0.25],	[0.20, 0.10]
				]
			},
			{
				"algo":	pollination.search,
				"params":[
					[0.7, 0.2]
				]
			},
			{
				"algo": fish.search,
				"params":[
					[0.125, 0.075],	[4, 2],		[0.125, 0.075],	[2, 1],
				]
			},
			{
				"algo": wolf.search,
				"params":[
					[0.5, 0.1],		[0.85, 0.15],	[0.5, 0.1],		[0.9, 0.075]	
				]
			}]
stage3 = [
			models.nelder_mead, models.powell, models.cg, models.bfgs, 
			models.lbfgsb, models.tnc, models.slsqp
			#models.cobyla, models.dogleg, models.trustncg
		]


# ---------- DONE, RUN SEARCH ALGORITHM -------------------------------------------------------------

nsga.set_global(globals())

pop = nsga.search(nsga_fn_evals, nsga_max_gens, nsga_pop_size, nsga_p_cross)
print("done!\n")


# ----------- VISUALIZATION -------------------------------------------------------------------------

pop.sort(key = lambda x: x["objectives"][0] + x["objectives"][1])

colors = {
	'bat':			'#020202',
	'bee':			'#bbaf22',
	'cuckoo':		'#ec994a',
	'firefly': 		'#ff0000',
	'fish':			'#3d7bff',
	'pollination':	'#ff3dff',
	'pso':			'#3dd6ff',
	'wolf':			'#e03e3e'
}

print('\033[4m' + "runtime         conv            algo    nm/ecm  gens    pop     params(algo-specific)" + '\033[0m')
for p in pop:
	r, r2, r3, params = nsga.decode(p["bitstring"])
	n = r2.__module__
	c = colors[n]

	print('\t'.join([str(round(x,8)) for x in p["objectives"]]).zfill(10), end="\t")
	print(f"{n}\t{r3.__name__}\t{params[2]}\t{len(params[3])}", end="\t")

	for i in params[4:]:
		print(round(i,3), end="\t")
	print()

	#plt.title("1+ep vs time")
	#plt.plot(p["objectives"][0],p["objectives"][1],c+"o")

'''
names = [ x[:-3] for x in os.listdir('algorithms') if x[-3:] == '.py']
prog = {}
for n in names:
	prog[n] = [x[n] for x in algs]
for n in names:
	plt.plot([l/nsga_pop_size for l in prog[n]],colors[n],label=n)
plt.title('algorithm counts vs generation')
plt.legend(loc='upper left')
'''


with open('output_populations.csv','w') as csvfile:
	model = models.models[sys.argv[1]]
	writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	writer.writerow([f"NSGA GENS: {nsga_max_gens}",f"NSGA POP: {nsga_pop_size}",f"NSGA CROSSOVER: {nsga_p_cross}",f"NSGA AVG EVALS: {nsga_fn_evals}",f"NSGA CROSSBREED: {nsga_cross_breed}"])
	writer.writerow([f"MODEL: {[x for x in models.models if models.models[x]==model][0]}", f"MODEL POP: {model_pop_count}", f"MODEL GENS: {model_generations}"])
	writer.writerow([])
	writer.writerow(["algorithm", "NM/ECM", "runtime", "AVG error", "RANDOM error", "best candidate's model parameters", "score of best"])
	for p in pop:
		r, r2, r3, params = nsga.decode(p["bitstring"])
		rep = list(params)
		sz = len(params[3])
		param_count = len(max(stage2, key = lambda x: len(x["params"]))["params"])
		pop_size = len(params[3])
		rep[3] = [						# generate population, pseudorandom
				[t[0] + t[1]*np.random.uniform(-1, 1) for t in model["estimates"]]
				for x in range(pop_size)]

		params = tuple(rep)
		newl = r2(*params)
		bst = min(newl, key = model['objective'])
		rs, con = r3(bst)
		sc = model['objective'](rs) / model['result']
		writer.writerow([r2.__module__, r3.__name__, p["objectives"][0],p["objectives"][1], sc, rs, model["objective"](rs)])

plt.show()