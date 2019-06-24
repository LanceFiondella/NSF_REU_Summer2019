from random import randrange, random

def onemax(bitstring):
	return len(bitstring.replace('0',''))	# count number of ones in bin string

def random_bitstring(num_bits):
	return str(bin(randrange(2**num_bits)))[2:].zfill(num_bits)	# generate num_bits length bin number

def binary_tournament(pop):
	popsize = len(pop)
	i = randrange(popsize)	# generate two random differing numbers
	j = (1 + i + randrange(popsize-1)) % popsize
	return pop[i] if (pop[i]["fitness"] > pop[j]["fitness"]) else pop[j]

def point_mutation(bitstring, rate = None):
	if rate == None:				# return a new bitstring with random bits flipped
		rate = 1 / len(bitstring)	# less bits are flipped if the string is longer
	t = "".join([(str(1-int(x)) if (random() < rate) else x ) for x in bitstring])
	return t

def crossover(p1, p2, rate):
	if random() > rate:				# get random mix of parent 1 and parent 2
		return p1

	point = randrange(1, len(p1)-1)
	return p1[:point] + p2[point:len(p1)]

def reproduce(selected, pop_size, p_cross, p_mutation):
	children = []
	for idx, p1 in enumerate(selected):
		p2 = selected[idx+1] if idx%2 == 0 else selected[idx-1]
		if p2 == selected[len(selected)-1]:
			p2 = selected[0]
		child = {"bitstring": crossover(p1["bitstring"], p2["bitstring"], p_cross)}
		child["bitstring"] = point_mutation(child["bitstring"], p_mutation)
		children.append(child)
		if len(children) >= pop_size:
			break
	return children

def search(max_gens, num_bits, pop_size, p_crossover, p_mutation):
	population = [{"bitstring": random_bitstring(num_bits)} for i in range(pop_size)]
	for c in population:
		c["fitness"] = onemax(c["bitstring"])
	population.sort(key = lambda x: x["fitness"])
	best = population[0]

	for gen in range(max_gens):
		selected = [binary_tournament(population) for x in range(pop_size)]
		children = reproduce(selected, pop_size, p_crossover, p_mutation)
		for c in children:
			c["fitness"] = onemax(c["bitstring"])
		population.sort(key = lambda x: x["fitness"])

		if children[0]["fitness"] >= best["fitness"]:
			best = children[0]
		population = children
		print(" > gen {}, best: {}, \t{}".format(gen, best["fitness"], best["bitstring"]))
		if best["fitness"] == num_bits:
			print("\tmaximized fitness in {} gens!".format(gen))
			break
	return best


num_bits = 64
max_gens = 100
pop_size = 100
p_crossover = 0.98
p_mutation = 1 / num_bits

best = search(max_gens, num_bits, pop_size, p_crossover, p_mutation)
print("done! best: f={}, s=\t{}".format(best["fitness"], best["bitstring"]))