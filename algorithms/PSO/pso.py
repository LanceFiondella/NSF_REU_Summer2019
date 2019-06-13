import numpy as np

def random_vector(dimensions):
	return [np.random.uniform(x[0],x[1]) for x in dimensions]

def search(objective, search_space, max_gen, population = None, pcount = 30, omega = 0.5, phi_p = 0.1, phi_g = 0.1):
	swarmbest = None								# have the swarm's best on hand
	pop = []

	if(population == None):							# generate random points if none given
		pop = [{"vector":random_vector(search_space)} for x in range(pcount)]
	else:
		pop = [{"vector":x} for x in population]	# convert an existing pop otherwise

	for p in pop:
		p["objective"] = objective(p["vector"])		# initialize objs / swarm best / velocities
		p["best_obj"] = p["objective"]
		p["best"] = p["vector"][:]

		if (swarmbest == None) or (p["best_obj"] < swarmbest["best_obj"]):
			swarmbest = p
		p["velocity"] = [np.random.uniform(-abs(x[1]-x[0]),abs(x[1]-x[0])) for x in search_space]

	for gen in range(max_gen):
		for p in pop:
			for d in range(len(search_space)):		# each dimension of each point, calc new velocity
				rp, rg = np.random.uniform(0, 1), np.random.uniform(0, 1)
				p["velocity"][d] =  omega * p["velocity"][d] \
									+ phi_p * rp * (p["best"][d] - p["vector"][d]) \
									+ phi_g * rg * (swarmbest["best"][d] - p["vector"][d])

			p["vector"] = list(np.add(p["vector"], p["velocity"]))	# move the point, find new obj
			p["objective"] = objective(p["vector"])
			if (p["objective"] < p["best_obj"]):				# if better than last, update it
				p["best"] = p["vector"][:]
				p["best_obj"] = p["objective"]
				if p["objective"] < swarmbest["best_obj"]:		# update swarm best if necessary
					swarmbest = p
	return [x["vector"] for x in pop]