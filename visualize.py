import importlib, sys

fname = sys.argv[1]
sys.argv[1] = sys.argv[2]	# to fix evalulate.py error

from evaluate import *
from matplotlib import pyplot as plt
from numpy import inf, linspace

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
	from matplotlib.ticker import MaxNLocator

	rtimecand4 = [0.07640242576599121, 0.07745361328125, 0.09600043296813965, 0.0966036319732666, 0.09758234024047852, 0.10510540008544922, 0.106201171875, 0.10864496231079102, 0.10931968688964844, 0.10937190055847168, 0.10947585105895996, 0.1095118522644043, 0.1095418930053711, 0.10954499244689941, 0.10956048965454102, 0.10956454277038574, 0.10959434509277344, 0.10968184471130371, 0.10969161987304688, 0.10971641540527344, 0.10988235473632812, 0.10996270179748535, 0.10997867584228516, 0.11004185676574707, 0.11006450653076172, 0.11011409759521484, 0.11020994186401367, 0.11023974418640137, 0.11032271385192871, 0.11033773422241211, 0.11035394668579102, 0.11045718193054199, 0.11066794395446777, 0.11076068878173828, 0.11079573631286621, 0.11084532737731934, 0.1110231876373291, 0.11107325553894043, 0.11120223999023438, 0.11128973960876465, 0.11170554161071777, 0.11206388473510742, 0.11228561401367188, 0.1125020980834961, 0.11357784271240234, 0.11498713493347168, 0.12031245231628418, 0.12042760848999023, 0.12056088447570801, 0.13097620010375977, 0.13276886940002441, 0.13333964347839355, 0.13547492027282715, 0.13582921028137207, 0.14239072799682617, 0.1459178924560547, 0.15554571151733398, 0.16669464111328125, 0.17372775077819824, 0.1743636131286621, 0.17766761779785156, 0.2211146354675293, 0.25800538063049316]
	acccand4 = [1.0903750354487216, 1.0903808275407032, 1.090386618820985, 1.0903898735458653, 1.0903971101377954, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904069348475065, 1.0904094610266384, 1.0904122511619623, 1.0904147808405815, 1.09041555368635, 1.0904226359266174, 1.0904790989515831, 1.0904814690453237, 1.0904890792098052, 1.090496511099383, 1.0904984651710723, 1.0905327766628612, 1.0905408447476979, 1.0905423383947757, 1.0905776332928316, 1.0905809486618745, 1.0906461542862251, 1.0906879915233707, 1.0907509926149948, 1.090776986104891, 1.0907979818626616, 1.0910889550308138, 1.0911653915007404, 1.0962200272386409]

	rtimecand5 = [0.07593417167663574, 0.07888269424438477, 0.08312654495239258, 0.08601021766662598, 0.0938112735748291, 0.09518599510192871, 0.09552001953125, 0.10059213638305664, 0.10292768478393555, 0.1068272590637207, 0.10757231712341309, 0.10897111892700195, 0.11124610900878906, 0.11167383193969727, 0.11457395553588867, 0.1177225112915039, 0.12216305732727051, 0.12415552139282227, 0.12457084655761719, 0.12666916847229004, 0.1272730827331543, 0.12879228591918945, 0.12937045097351074, 0.13037467002868652, 0.1339709758758545, 0.13399124145507812, 0.13776230812072754, 0.13939809799194336, 0.1409015655517578, 0.14179229736328125, 0.14284086227416992, 0.14334940910339355, 0.14344573020935059, 0.14523673057556152, 0.1475837230682373, 0.14870882034301758, 0.14881134033203125, 0.15183758735656738, 0.1524193286895752, 0.15439701080322266, 0.15474700927734375, 0.15784096717834473, 0.16099023818969727, 0.1625685691833496, 0.16866731643676758, 0.1694014072418213, 0.17367172241210938, 0.1778569221496582, 0.17968392372131348, 0.17981982231140137, 0.18387055397033691, 0.187483549118042, 0.18751144409179688, 0.20161747932434082, 0.2030773162841797, 0.2049083709716797, 0.22011446952819824, 0.24331235885620117, 0.354567289352417, 0.3706674575805664, 0.37616467475891113, 0.39525294303894043, 0.4162893295288086]
	acccand5 = [1.0903790906111508, 1.0903819936195023, 1.0903824908478996, 1.0903901844696193, 1.0903947759707846, 1.0903957788398675, 1.0903999835233105, 1.0904002228248015, 1.0904116145658567, 1.0904216147679313, 1.0904226615107933, 1.0904308740651507, 1.0904334244346883, 1.0904458169146987, 1.0904503380745276, 1.0904614292613874, 1.0904703852040183, 1.0904743750411123, 1.0904766915745716, 1.0904813061685588, 1.0904849940537613, 1.0904916463837016, 1.0904933949670432, 1.090494276172457, 1.090494408680382, 1.0905377546102293, 1.0905465665757714, 1.0905567118264174, 1.0905745516909655, 1.0906032949437816, 1.0906052767142784, 1.0906472123999253, 1.0906491407269596, 1.0906600763291914, 1.0906631182523019, 1.0906945298233155, 1.0907185732288114, 1.0907548346620883, 1.090782947276038, 1.0907871276690895, 1.0908798790190748, 1.0909034318992528, 1.0909069859998728, 1.091170869544874, 1.091918278823773, 1.098151928998124, 1.1567527001618703, 1.1714063020332055, 1.2011582246912933, 1.2251065980510578, 1.228897848591283, 1.2389077966316282, 1.2578746068641604, 1.3803705750586748, 1.47296597893819, 1.5648255692913853, 1.6115030818093536, 2.2287581489416204]#, inf, inf, inf, inf, inf]
	
	acccand4 = [t/min(acccand4) for t in acccand4]	# fix wrong MLE issue
	acccand5 = [t/min(acccand5) for t in acccand5]	# fix wrong MLE issue

	hcount = 30

	ax = plt.figure().gca()
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	plt.hist(rtimecand4,hcount, color='#000000',label='Design 4')
	plt.hist(rtimecand5, hcount, color='#AAAAAA88',label='Design 5')
	plt.ylabel('Count')
	plt.xlabel('Runtime (s)')
	plt.legend(loc='best')
	plt.show()

	hcount = linspace(1, 2.5, hcount)

	ax = plt.figure().gca()
	plt.hist(acccand4,hcount, color='#000000',label='Design 4')
	plt.hist(acccand5, hcount, color='#AAAAAA88',label='Design 5')
	plt.ylabel('Count')
	plt.xlabel('Accuracy ($1+\epsilon$)')
	plt.legend(loc='best')
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





