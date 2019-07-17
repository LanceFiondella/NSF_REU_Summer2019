import numpy as np

def random_vector(dimensions):
	return np.array([np.random.uniform(x[0],x[1]) for x in dimensions])

def search(objective, search_space, max_gen, population = None, omega = 0.5, phi_p = 0.1, phi_g = 0.1):
	swarmbest = None								# have the swarm's best on hand
	pop = []

	if(population == None):							# generate random points if none given
		pop = [{"vector":random_vector(search_space)} for x in range(pcount)]
	else:
		pop = [{"vector":np.array(x)} for x in population]	# convert an existing pop otherwise

	for p in pop:
		p["objective"] = objective(p["vector"])		# initialize objs / swarm best / velocities
		p["best_obj"] = p["objective"]
		p["best"] = p["vector"].copy()

		if (swarmbest == None) or (p["best_obj"] < swarmbest["best_obj"]):
			swarmbest = p
		p["velocity"] = np.array([np.random.uniform(-abs(x[1]-x[0]),abs(x[1]-x[0])) for x in search_space])

	for gen in range(max_gen):
		for p in pop:
			p['velocity'] = omega * p['velocity'] + \
							phi_p * (p["best"] - p["vector"]) * np.array([np.random.uniform() for x in p['vector']]) + \
							phi_g * (swarmbest["best"] - p["vector"]) * np.array([np.random.uniform() for x in p['vector']])

			p["vector"] = p['vector'] + p['velocity']
			
			for i in range(len(p["vector"])):		# keep it inside search space
				p["vector"][i] = min(p["vector"][i], search_space[i][1])
				p["vector"][i] = max(p["vector"][i], search_space[i][0])
				
			p["objective"] = objective(p["vector"])

			if (p["objective"] < p["best_obj"]):				# if better than last, update it
				p["best"] = p["vector"].copy()
				p["best_obj"] = p["objective"]
				if p["objective"] < swarmbest["best_obj"]:		# update swarm best if necessary
					swarmbest = p
	return [x["vector"] for x in pop]

if __name__ == "__main__":
	obj = lambda x: x[0]**2 + x[1]**2
	ss = [[-3,3],[-3,3]]
	pop = [[np.random.uniform(),np.random.uniform()] for i in range(25)]
	print(search(obj, ss, 30, pop))