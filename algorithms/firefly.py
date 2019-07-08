import numpy as np

def attractiveness(a, b, absorption):		# relative attractiveness between two fireflies
	radsq = sum([ (a[i] - b[i])**2 for i in range(len(a))])	# dist squared from a<->b
	return 1/(1 + radsq * absorption)		# assume B0 = 1, yang pg 113


def random_vector(search_space):
	return np.array([np.random.uniform(x[0], x[1]) for x in search_space])


def gaussian_random(scalar, search_space, generation):
	alpha = scalar ** generation
	return np.array([(np.random.random_sample() - 0.5) * alpha for x in range(len(search_space))])


def search(objective, search_space, max_gens, pop = None, pop_count = 30, absorption = 1, randomness = 0.99):

	if pop == None:
		pop = [{"vector":random_vector(search_space)} for x in range(pop_count)]
	else:
		pop = [{"vector":np.array(x)} for x in pop]

	for p in pop:
		p["objective"] = objective(p["vector"]) #intensity(p["vector"], absorption)

	last = 0
	for gen in range(max_gens):
		for p in pop:
			for n in pop:
				if n["objective"] <= p["objective"]:
										# calculate random movement, decreases with time	(8.13)
										# attractiveness between two fireflies 				(8.4)
					attr = attractiveness(p["vector"], n["vector"], absorption)
										# calculate new position for firefly 				(8.11)
					p["vector"] = 	p["vector"] + \
									gaussian_random(randomness, search_space, gen) + \
									attr * (n["vector"] - p["vector"])

					for i in range(len(p["vector"])):		# keep it inside search space
						p["vector"][i] = min(p["vector"][i], search_space[i][1])
						p["vector"][i] = max(p["vector"][i], search_space[i][0])

					p["objective"] = objective(p["vector"])	# re-calculate objective value
					
	return [x["vector"] for x in pop] 		# return the converged population

if __name__ == "__main__":
	obj = lambda x: x[0]**2 + x[1]**2
	search_space = [[0,1],[0,1]]
	max_gens = 20
	pop = [[0.1,0] for x in range(30)]
	print(search(obj, search_space, max_gens, pop))