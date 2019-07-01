import numpy as np

def search(objective, search_space, max_gen, pop = None, pop_count = 30, look_pct = 0.33, exp_pct = 0.33, wb = 0.5, wg = 0.5, rmult = 0.2):

	pop = [{"vector": np.array(x), 'best_objective':float('inf')} for x in pop]
	
	for gen in range(max_gen):

		for p in pop:								# fit the bees for sorting, set up individual best if necessary
			p["objective"] = objective(p["vector"])	# (will always update best on 1st iteration)
			if (gen == 0) or (p["objective"] < p["best_objective"]):
				p["best_objective"] = p["objective"]
				p["best_vector"] = p["vector"].copy()

		pop.sort(key = lambda x: x["objective"])

		best = pop[0]["vector"].copy()

		for i, p in enumerate(pop):

			if i < exp_pct*len(pop):				# "experienced" bees
													# move in average to best and local best
				p["vector"] = 	p["vector"] + \
								wg * np.random.uniform() * (best - p["vector"]) + \
								wb * np.random.uniform() * (p["best_vector"] - p["vector"])

			elif i < (look_pct+exp_pct)*len(pop):	# "onlooker" bees

				elite = np.random.choice(pop[:int(len(pop)*exp_pct)])
													# get a random experienced bee
				p["vector"] = 	p["vector"] + \
								wg * np.random.uniform() * (elite["vector"] - p["vector"])
													# move it towards elite bee

			else:									# "scout" bees
				nsearch = [np.random.uniform(x[0],x[1])*rmult*np.random.uniform() for x in search_space]
													# get pseudo-radius walk in search space
				p["vector"] = 	p["vector"] + np.array(nsearch)
		if __name__ == "__main__":
			print(f" > iteration {gen+1}")
	return [x["vector"] for x in pop]

if __name__ == "__main__":
	objective = lambda x: x[0]**2 + x[1]**2
	search_space = [[-5,5],[-5,5]]
	max_gen = 10
	pop = [[np.random.uniform(-5,5),np.random.uniform(-5,5)] for x in range(75)]
	npop = search(objective, search_space, max_gen, pop)
	print("done!")
	print(min(npop, key = lambda x: objective(x)))