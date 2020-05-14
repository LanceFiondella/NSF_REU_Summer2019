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

if False:	
	from matplotlib.ticker import MaxNLocator

	rtimecand5 = [0.2709615230560303, 0.29187822341918945, 0.30053043365478516, 0.32530927658081055, 0.3485677242279053, 0.3488740921020508, 0.35234522819519043, 0.35310959815979004, 0.3610849380493164, 0.3618941307067871, 0.3670201301574707, 0.37003540992736816, 0.37067437171936035, 0.37110352516174316, 0.3758275508880615, 0.3796424865722656, 0.3940620422363281, 0.3956789970397949, 0.3985741138458252, 0.40567660331726074, 0.40799641609191895, 0.40840768814086914, 0.4084334373474121, 0.4120018482208252, 0.41673731803894043, 0.4201350212097168, 0.4228050708770752, 0.4301881790161133, 0.431443452835083, 0.43282556533813477, 0.4329864978790283, 0.43717432022094727, 0.4372079372406006, 0.4440462589263916, 0.44664621353149414, 0.44736266136169434, 0.45020103454589844, 0.4566836357116699, 0.4658060073852539, 0.4670283794403076, 0.46713829040527344, 0.472642183303833, 0.48274850845336914, 0.49042844772338867, 0.4980041980743408, 0.509042501449585, 0.51418137550354, 0.5406920909881592, 0.5436422824859619, 0.5446023941040039, 0.5482847690582275, 0.564993143081665, 0.5697484016418457, 0.5867340564727783, 0.5899932384490967, 0.6013312339782715, 0.629901647567749, 0.6440277099609375, 0.6571774482727051, 0.7277491092681885, 0.7368276119232178, 0.8336653709411621, 0.9046330451965332]
	acccand5 = [1.0896114748465147, 1.0901623464659542, 1.090278574965895, 1.0902906176745946, 1.0903019984037052, 1.0903141508804424, 1.0903230601626184, 1.090344521456849, 1.0903543664866906, 1.0903574509098177, 1.0903604255828054, 1.0903615954664843, 1.0903623143048502, 1.0903665621580267, 1.090367477330004, 1.0903680519909544, 1.0903689469739082, 1.0903693732467747, 1.09036938203206, 1.0903694705896132, 1.0903694841375775, 1.0903694996591735, 1.0903695495413313, 1.0903696172520874, 1.0903696324817465, 1.0903696573595647, 1.0903696770516893, 1.0903696836228562, 1.0903697206704144, 1.0903697308864764, 1.0903697328819444, 1.0903697486998456, 1.0903697506787546, 1.0903697538009165, 1.0903697605788345, 1.090369761928695, 1.0903697628920437, 1.0903697638136864, 1.090369766322474, 1.090369766862013, 1.0903697669356704, 1.0903697676118904, 1.0903697707407578, 1.0903697708448548, 1.090369772251216, 1.0903697733664752, 1.0903697742971166, 1.090369774881361, 1.0903697756058195, 1.0903697757287136, 1.0903697771411727, 1.090369777846641, 1.0903697779019312, 1.0903697803279935, 1.0903697803437766, 1.0903697804649364, 1.0903697853663028, 1.0903697862133246, 1.0903697862368806, 1.0903697948913404, 1.0903698042925645, 1.0903704984772868, 1.0903720436572437]

	rtimecand6 = [0.08158516883850098, 0.08311724662780762, 0.08816695213317871, 0.0912632942199707, 0.09194636344909668, 0.0966944694519043, 0.09827780723571777, 0.0991816520690918, 0.09935569763183594, 0.1014862060546875, 0.10268521308898926, 0.1039278507232666, 0.10414528846740723, 0.10616636276245117, 0.10883450508117676, 0.10922837257385254, 0.11136221885681152, 0.11213397979736328, 0.11227560043334961, 0.11514830589294434, 0.11516427993774414, 0.11619138717651367, 0.12131166458129883, 0.12423372268676758, 0.12533330917358398, 0.12653851509094238, 0.12966275215148926, 0.13474822044372559, 0.13816452026367188, 0.13913369178771973, 0.13977265357971191, 0.14589381217956543, 0.1470644474029541, 0.148284912109375, 0.14977097511291504, 0.1520223617553711, 0.1531527042388916, 0.15474987030029297, 0.15591096878051758, 0.15665268898010254, 0.16195011138916016, 0.16287994384765625, 0.16309452056884766, 0.17009568214416504, 0.17123961448669434, 0.17483806610107422, 0.17746424674987793, 0.17970943450927734, 0.18782353401184082, 0.188673734664917, 0.18871259689331055, 0.18936729431152344, 0.19984817504882812, 0.20831680297851562, 0.21533513069152832, 0.22229433059692383, 0.23248553276062012, 0.24061036109924316, 0.24249529838562012, 0.25199246406555176, 0.3464624881744385, 0.3719172477722168, 0.3816647529602051]
	acccand6 = [1.0903783471083341, 1.0903791533252294, 1.0903852431880316, 1.0903939245974983, 1.0903960371727646, 1.0903985457382368, 1.0904031574683624, 1.0904102331451762, 1.0904117129752415, 1.0904135810881441, 1.090427789280318, 1.0904367546454015, 1.0904403225108668, 1.0904480029392463, 1.0904486518724792, 1.090453367637103, 1.0904559584961477, 1.0904652026526347, 1.090475661381195, 1.0904763959355541, 1.0904961516832885, 1.0905006536182975, 1.0905378952675417, 1.090548817685427, 1.0905576845904705, 1.0905655596731734, 1.0906021898804965, 1.0906325136477446, 1.0906340877365555, 1.0906356915678457, 1.0906770650949469, 1.090715161918913, 1.0907653134292836, 1.0907920173954675, 1.0917117913329841, 1.0918245471611157, 1.0919344618848457, 1.09219935767643, 1.09353570355204, 1.095935145735172, 1.0992229646580207, 1.1063702953129158, 1.1107643290553257, 1.1300339255447425, 1.1585065626908453, 1.166935849011623, 1.1845229911713853, 1.1867847397745699, 1.196960908617409, 1.2087156309943283, 1.2095791384261214, 1.2351471291309977, 1.2370977071661131, 1.2391076263492464, 1.268228432555127, 1.3105179653129853, 1.3680040688495407, 1.4184294260046197, 1.7470064187186969, 2.404843673055993] #,inf,inf,inf]

	hcount = 30

	ax = plt.figure().gca()
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	plt.hist(rtimecand5,hcount, color='#404040')
	plt.ylabel('Number of runs')
	plt.xlabel('Runtime (s)')
	plt.show()

	ax = plt.figure().gca()
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	plt.ticklabel_format(useOffset=False)
	plt.hist(acccand5, hcount, color='#404040')
	plt.ylabel('Number of runs')
	plt.xlabel('Accuracy ($1+\epsilon$)')
	plt.show()

	ax = plt.figure().gca()
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	plt.hist(rtimecand6, hcount, color='#404040')
	plt.ylabel('Number of runs')
	plt.xlabel('Runtime (s)')
	plt.show()

	ax = plt.figure().gca()
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	plt.hist(acccand6, hcount, color='#404040')
	plt.ylabel('Number of runs')
	plt.xlabel('Accuracy ($1+\epsilon$)')
	plt.show()


if True:	# print all candidates
	plt.ticklabel_format(useOffset=False)
	for p in pops[-1]:
		a,b,c,d = decode(p['bitstring'])
		p = p['objectives']
		#print(a)
		#if b == None:
		plt.scatter(p[0], p[1], color='k' if b != None else '#808080',marker='.')

	plt.xlim([0, .225])#plt.xlim([0 ,0.225])
	plt.ylim([1-0.1, 4.25])#plt.ylim([1-0.025, 3.75])
	plt.xlabel('Runtime (s)',fontsize=12)
	plt.ylabel('Accuracy ($1+\epsilon$)',fontsize=12)
	#plt.title('covariate dataset 2')
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





