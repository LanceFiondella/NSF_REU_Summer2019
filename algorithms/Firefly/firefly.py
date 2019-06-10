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

def search(max_gens, search_space, pop_count, absorption, randomness):
	pop = [{"vector":random_vector(search_space)} for x in range(pop_count)]
	for p in pop:
		p["intensity"] = intensity(p["vector"], absorption)

	for gen in range(max_gens):
		for i in range(len(pop)):
			for j in range(i):
				if pop[j]["intensity"] > pop[i]["intensity"]:
										# calculate random movement, decreases with time	(8.13)
					rands = gaussian_random(randomness, search_space, gen)
										# attractiveness between two fireflies 				(8.4)
					attr = attractiveness(pop[i]["vector"], pop[j]["vector"], absorption)
										# calculate new position for firefly 				(8.11)
					term = np.multiply(attr, np.subtract(pop[j]["vector"], pop[i]["vector"]))
					term = np.add(term, rands)

					pop[i]["vector"] = np.add(pop[i]["vector"], term)
					pop[i]["intensity"] = intensity(pop[i]["vector"], absorption)

		best = max(pop, key = lambda x: x["intensity"])
		print(" > gen {}: best {}, {}".format(1+gen, best["vector"], best["intensity"]))
	return best

problem_size = 1
search_space = [[-10, 10] for x in range(problem_size)]
pop_count = 25
max_gens = 25
absorption = 1	# see pg 115, 0.001 to 1000
randomness = 0.95	# 0.95 -> 0.99

best = search(max_gens, search_space, pop_count, absorption, randomness)
print("done!")