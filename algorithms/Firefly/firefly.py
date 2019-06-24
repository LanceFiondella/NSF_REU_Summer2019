import numpy as np

def intensity(vector, absorption):			# problem to maximize
	radsq = sum([x**2 for x in vector])		# r = sqrt(x^2 + y^2): e^-yr^2
	return 1/np.exp(radsq * absorption)		# see yang pg 113 - assume original intensity at src is 1

def attractiveness(a, b, absorption):		# relative attractiveness between two fireflies
	radsq = sum([ (a[i] - b[i])**2 for i in range(len(a))])	# dist squared from a<->b
	return 1/(1 + radsq * absorption)		# assume B0 = 1, yang pg 113


def random_vector(search_space):
	return [np.random.uniform(x[0], x[1]) for x in search_space]


def gaussian_random(scalar, search_space, generation):
	alpha = scalar ** generation
	return [(np.random.random_sample() - 0.5) * alpha for x in range(len(search_space))]


def search(objective, search_space, max_gens, pop = None, pop_count = 30, absorption = 1, randomness = 0.99):

	if pop == None:
		pop = [{"vector":random_vector(search_space)} for x in range(pop_count)]
	else:
		pop = [{"vector":x} for x in pop]

	for p in pop:
		p["objective"] = objective(p["vector"]) #intensity(p["vector"], absorption)

	last = 0
	for gen in range(max_gens):
		for p in pop:
			for n in pop:
				if n["objective"] <= p["objective"]:
										# calculate random movement, decreases with time	(8.13)
					rands = gaussian_random(randomness, search_space, gen)
										# attractiveness between two fireflies 				(8.4)
					attr = attractiveness(p["vector"], n["vector"], absorption)
										# calculate new position for firefly 				(8.11)
					term = np.multiply(attr, np.subtract(n["vector"], p["vector"]))
					term = np.add(term, rands)

					p["vector"] = np.add(p["vector"], term)
					for i in range(len(p["vector"])):
						if p["vector"][i] < search_space[i][0]:
							p["vector"][i] = search_space[i][0]
						elif p["vector"][i] > search_space[i][1]:
							p["vector"][i] = search_space[i][1]

					p["objective"] = objective(p["vector"]) #intensity(p["vector"], absorption)
					if np.isnan(p["objective"]):
						p["objective"] = float('inf')
	return [x["vector"] for x in pop] 		# return the converged population

if __name__ == "__main__":
	obj = lambda x: x[0]**2 + x[1]**2
	search_space = [[0,1],[0,1]]
	max_gens = 20
	pop = [[0.1,0] for x in range(30)]
	print(search(obj, search_space, max_gens, pop))