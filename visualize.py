import importlib, sys

fname = sys.argv[1]
sys.argv[1] = sys.argv[2]	# to fix evalulate.py error

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
for p in pop:
	r, r2, r3, params = decode(p["bitstring"])
	n = r2.__module__ if r2 != None else "NONE"
	c = colors[n]

	print(sep.join([str(round(x,12)).ljust(12, ' ') for x in p["objectives"]]), end=sep)

	print(f"{n[:6] if r2 != None else 'NONE'}{sep}{r3.__name__[:6] }{sep}{params[2] if r2 != None else 'n/a'}{sep}{len(params[3])if r2 != None else 'n/a'}", end=sep)
	if r2 != None:
		for i in params[4:]:
			print(round(i,3), end=sep)
	print()

	#print(p["bitstring"], end=sep)
	#print(r,p['objectives'][0], p['objectives'][1])
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

if True:	
	rtimecand5 = [0.41565752029418945, 0.42770886421203613, 0.48233723640441895, 0.48890209197998047, 0.4921989440917969, 0.4974536895751953, 0.5026638507843018, 0.5032916069030762, 0.503307580947876, 0.5107781887054443, 0.5132715702056885, 0.5152654647827148, 0.5171256065368652, 0.5287964344024658, 0.5554893016815186, 0.558380126953125, 0.5662267208099365, 0.5746846199035645, 0.5778200626373291, 0.5778775215148926, 0.5859389305114746, 0.5920743942260742, 0.593203067779541, 0.5975706577301025, 0.6036665439605713, 0.6037688255310059, 0.6043570041656494, 0.6045722961425781, 0.6099808216094971, 0.6136839389801025, 0.6148571968078613, 0.6152803897857666, 0.6244308948516846, 0.638129711151123, 0.6408987045288086, 0.6476624011993408, 0.6516854763031006, 0.6528236865997314, 0.6574647426605225, 0.6631057262420654, 0.6653048992156982, 0.6772494316101074, 0.6818523406982422, 0.68996262550354, 0.6909756660461426, 0.6918315887451172, 0.6921651363372803, 0.6935269832611084, 0.6954498291015625, 0.695458173751831, 0.6995861530303955, 0.7117934226989746, 0.7154228687286377, 0.7182271480560303, 0.734447717666626, 0.7448215484619141, 0.7706923484802246, 0.7763869762420654, 0.8113677501678467, 0.8477857112884521, 0.8850581645965576, 0.8874838352203369, 0.9971272945404053]
	acccand5 = [1.086908748160055, 1.087042830473172, 1.0903056921915242, 1.090346277578079, 1.0903548446496054, 1.0903549971931217, 1.0903626432484625, 1.0903648462258329, 1.0903648999493383, 1.0903661454343714, 1.090367141403077, 1.0903671765538925, 1.0903679827897828, 1.09036816978423, 1.090368267189484, 1.0903687506702322, 1.0903688313012878, 1.0903691981230743, 1.0903693541163715, 1.0903694682654488, 1.0903695341929038, 1.0903695693740652, 1.0903696476738456, 1.0903696792709365, 1.0903696828161225, 1.0903697135408121, 1.090369719869289, 1.0903697206721241, 1.0903697244177113, 1.0903697260138099, 1.0903697370510548, 1.0903697394563037, 1.0903697593325128, 1.090369762996368, 1.0903697634078637, 1.0903697644725472, 1.0903697696379657, 1.09036977014077, 1.0903697713397278, 1.0903697714797986, 1.090369772742378, 1.0903697728341213, 1.090369773169276, 1.0903697732204303, 1.0903697742444554, 1.090369774547387, 1.0903697748141075, 1.0903697758978501, 1.0903697762196287, 1.0903697766653588, 1.0903697776398902, 1.09036977869939, 1.0903697791462665, 1.0903697795903413, 1.0903697801632153, 1.0903697815164854, 1.090369784018167, 1.0903697944855175, 1.0903697984584813, 1.0903699047834974, 1.0903701662905407, 1.090370631903534, 1.0903709959527574]

	rtimecand6 = [0.10234451293945312, 0.1085364818572998, 0.1240546703338623, 0.12560439109802246, 0.1271359920501709, 0.12818527221679688, 0.13466691970825195, 0.13474631309509277, 0.14086031913757324, 0.1459980010986328, 0.16624927520751953, 0.16974878311157227, 0.17142367362976074, 0.17654633522033691, 0.18280553817749023, 0.18461871147155762, 0.18686151504516602, 0.1877453327178955, 0.18815183639526367, 0.19000506401062012, 0.19470500946044922, 0.1961374282836914, 0.20194745063781738, 0.20200514793395996, 0.20848846435546875, 0.21070408821105957, 0.21305346488952637, 0.21350431442260742, 0.2142937183380127, 0.2153167724609375, 0.21659135818481445, 0.22159409523010254, 0.22491455078125, 0.2250349521636963, 0.23126220703125, 0.23153471946716309, 0.2317497730255127, 0.23333311080932617, 0.238175630569458, 0.24240589141845703, 0.24341797828674316, 0.24657249450683594, 0.24840617179870605, 0.24991250038146973, 0.2533450126647949, 0.26102781295776367, 0.2626686096191406, 0.2663140296936035, 0.2736358642578125, 0.277071475982666, 0.27840471267700195, 0.2814974784851074, 0.28175950050354004, 0.2852799892425537, 0.2929389476776123, 0.2973806858062744, 0.316784143447876, 0.32663893699645996, 0.32793545722961426, 0.34613871574401855, 0.3620922565460205, 0.430128812789917, 0.5996921062469482]
	acccand6 = [1.0903818862766677, 1.0903968177398562, 1.0903984640582611, 1.0903999822354604, 1.0904052403900957, 1.0904075483415892, 1.0904096169441828, 1.0904231275933, 1.0904275093118274, 1.090434421422838, 1.0904419193905959, 1.090449672837005, 1.0904671506685717, 1.090482789066326, 1.0904944844830318, 1.0905058193910204, 1.090511833163876, 1.0905210350708092, 1.0905276712650849, 1.0905323436141963, 1.090555408308006, 1.0905673986750914, 1.0905819544900361, 1.09060199594318, 1.0906180002418095, 1.090633279652438, 1.0906455894640756, 1.0906499778152272, 1.0907001593816628, 1.0907009111601904, 1.0907429709804128, 1.0907717152365985, 1.0907755703786721, 1.0909272280943367, 1.0910379997129822, 1.0910760769903203, 1.0911178119552676, 1.091194013826477, 1.0937231803505545, 1.0950081001727539, 1.09738234939417, 1.0993441198506217, 1.1021581669263267, 1.1411420959650445, 1.1607939400180565, 1.1827231740742747, 1.198741129315511, 1.2181982131371962, 1.2357182109940463, 1.239904092558028, 1.264048735542887, 1.276112643322295, 1.369152220967698, 1.413872652935625, 1.4244498463020978, 1.4434016378840777, 1.59648304586574, 1.5972431850766957, 1.6117962337148284, 1.631405458999721, 2.512954812202093]

	plt.hist(rtimecand5)
	plt.ylabel('Number of runs')
	plt.xlabel('Runtime')
	plt.show()

	plt.hist(acccand5)
	plt.ylabel('Number of runs')
	plt.xlabel('Accuracy (1+\epsilon)')
	plt.show()

	plt.hist(rtimecand6)
	plt.ylabel('Number of runs')
	plt.xlabel('Runtime')
	plt.show()

	plt.hist(acccand6)
	plt.ylabel('Number of runs')
	plt.xlabel('Accuracy (1+\epsilon)')
	plt.show()


if True:	# print all candidates
	plt.ticklabel_format(useOffset=False)
	for p in pops[-1]:
		a,b,c,d = decode(p['bitstring'])
		p = p['objectives']
		#print(a)
		#if b == None:
		plt.scatter(p[0], p[1], color='k' if b != None else 'b',marker='.')

	plt.xlim([0, .25])#plt.xlim([0 ,0.225])
	plt.ylim([1-0.025, 5])#plt.ylim([1-0.025, 3.75])
	plt.xlabel('Runtime (s)',fontsize=12)
	plt.ylabel('Accuracy ($1+\epsilon$)',fontsize=12)
	plt.title('covariate dataset 2')
	plt.show()

else:
	#print weibull data
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
			set_ind = int(ind/3)
			dset = dsets[set_ind]
			if dset == "SYS1" and False:	# option to skip sys1
				continue
			if ind%3 == 2 and True:		# option to skip 2nd bfgs
				continue
			#print(f"{dset} {'bee' if ind%3 == 0 else 'bfgs'} {p['objectives'][0]} {p['objectives'][1]}")
			color = colors[ind%3]
			ecolor = None #if color == 'k' else 'k'
			plt.scatter(p['objectives'][0], p['objectives'][1], marker=marks[set_ind], color=color, edgecolors=ecolor, label=dsets[set_ind] if ind%3==0 else None)
			plt.ylim([1-0.000125, 1.002125])
			
	#plt.xlim([0,0.015])
	plt.xlabel('Runtime (s)',fontsize=12)
	plt.ylabel('Accuracy ($1+\epsilon$)',fontsize=12)

	plt.legend(loc='best',scatterpoints = 1)
	plt.show()





