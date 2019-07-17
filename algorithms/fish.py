import numpy as np

def objective_function(u):
	"""Objective function used for testing. 3D sphere function.

	Keyword arguments:
	u (list) -- (x, y) values

	Returns:
	float -- value of objective function evaluated at (x, y)
	"""
	#z = (1 - u[0])**2 + 100 * (u[1] - u[0]**2)**2 + (1 - u[2])**2   # 3D Rosenbrock function
	z = u[0]**2 + u[1]**2
	return z

def calc_visual_scope(current_fish, visual, problem_size, search_space):
	"""Calculates the visual scope bounds of a specified fish. Bounds are a simple
	minimum and maximum for each dimesion of problem, giving the visual a square shape.

	Keyword arguments:
	current_fish (dict) -- fish whose visual scope is calculated
	visual (float) -- distance from fish to edge of visual (radius, if circular)
	problem_size (int) -- dimension of problem (number of independent variables)
	search_space (list) -- search space of problem, bounds for each dimension
	"""
	# (x - center_x)^2 + (y - center_y)^2 < radius^2
	pos = current_fish["position"]
	# this creates the visual
	current_fish["visual"] = np.array([[pos[j] - visual, pos[j] + visual] for j in range(problem_size)])
	for j in range(problem_size):
		if (current_fish["visual"][j][0] < search_space[j][0]):
			current_fish["visual"][j][0] = search_space[j][0]
		elif (current_fish["visual"][j][1] > search_space[j][1]):
			current_fish["visual"][j][1] = search_space[j][1]
	# !!! check for bounds here or when moves made???

def is_in_visual(fish_i, fish_j, problem_size):
	"""Determines if fish_j is in the visual range of fish_i.

	Keyword arguments:
	fish_i (dict) -- fish whose visual range will be used
	fish_j (dict) -- fish checked to see if in fish_i's visual range
	problem_size (int) -- dimension of problem (number of independent variables)

	Returns:
	bool -- True if in visual, False if not in visual
	"""
	#i_pos = fishi["position"]
	vis = fish_i["visual"]
	j_pos = fish_j["position"]
	if ((j_pos[k] > vis[k][0] and j_pos[k] < vis[k][1]) for k in range(problem_size)):
		return True
	else:
		return False

def which_fish_in_visual(pop, i, problem_size):
	'''Determines which fish are in the visual of a given fish.

	Keyword arguments:
	pop (list) -- entire fish population, elements are dictionaries representing individual fish
	i (int) -- current number of main loop iteration
	problem_size (int) -- dimension of problem (number of independent variables)

	Returns:
	list -- indexes of fish in the visual
	'''
	fish_in_vis = []   # number of fish in visual, list of indexes of fish in visual
	for j in range(len(pop)):
		j_in_visual = is_in_visual(pop[i], pop[j], problem_size)
		if (j_in_visual and i != j):
			fish_in_vis.append(j)    # keep track of which fish are in visual, used to calc centroid
	return fish_in_vis

def calc_centroid(pop, fish_in_vis, num, problem_size):
	'''Calculates the centroid (average position) of fish in the visual of a given fish.

	Keyword arguments:
	pop (list) -- entire fish population, elements are dictionaries representing individual fish
	fish_in_vis (list) -- current number of main loop iteration
	problem_size (int) -- dimension of problem (number of independent variables)

	Returns:
	list -- indexes of fish in the visual
	'''
	# returns centroid position
	pos_sum = np.zeros(problem_size)
	for j in range(problem_size):
		pos_sum[j] = np.sum(pop[index]["position"][j] for index in fish_in_vis)
	return np.array(pos_sum / num)

def move_state(current_fish, search_space, problem_size, objective):
	'''Executes the move state behavior. Fish moves to a random position in search
	space if it has better fitness.

	Keyword arguments:
	current_fish (dict) -- fish whose visual scope is calculated
	search_space (list) -- search space of problem, bounds for each dimension
	problem_size (int) -- dimension of problem (number of independent variables)
	objective (function) - objective function used to evaluate population fitness
	'''
	# move to random position in search space if it has better fitness
	rand_pos = np.zeros(problem_size)
	for j in range(problem_size):
		rand_pos[j] = np.random.uniform(search_space[j][0], search_space[j][1])
	if (objective(rand_pos) < current_fish["fitness"]):
		current_fish["position"] = rand_pos
		current_fish["fitness"] = objective(current_fish["position"])

def prey_state(pop, current_fish, fish_in_vis, try_number, step, objective):
	'''Executes the prey state behavior. Fish moves toward random fish in visual.

	Keyword arguments:
	pop (list) -- entire fish population, elements are dictionaries representing individual fish
	current_fish (dict) -- fish whose visual scope is calculated
	fish_in_vis (list) -- current number of main loop iteration
	try_number (int) -- number of times prey behavior is executed
	step (float) -- max step size fish can move
	objective (function) - objective function used to evaluate population fitness
	'''
	rand_fish = {}
	for i in range(try_number):
		rand_index = np.random.choice(fish_in_vis)
		rand_in_vis = pop[rand_index]
		rand_fish["position"] = np.add(current_fish["position"], rand_in_vis["position"])
		rand_fish["fitness"] = objective(rand_fish["position"])
		if (rand_fish["fitness"] < current_fish["fitness"]):
			rand_minus_current = np.subtract(rand_fish["position"], current_fish["position"])
			magnitude = np.absolute(np.subtract(rand_fish["position"], current_fish["position"]))
			direction = np.divide(rand_minus_current, magnitude)
			current_fish["position"] = np.add(current_fish["position"], np.random.random() * step * direction)
		else:
			current_fish["position"] = np.add(current_fish["position"], np.random.uniform(-1.0, 1.0) * step)
	current_fish["fitness"] = objective(current_fish["position"])

def swarm_state(current_fish, center_pos, step, objective):
	'''Executes the swarm state behavior. Fish moves toward cener of fish in visual.

	Keyword arguments:
	current_fish (dict) -- fish whose visual scope is calculated
	center_pos (list) -- position of center of fish in visual
	step (float) -- max step size fish can move
	objective (function) - objective function used to evaluate population fitness
	'''
	center_minus_current = np.subtract(center_pos, current_fish["position"])
	magnitude = np.absolute(np.subtract(center_pos, current_fish["position"]))
	direction = np.divide(center_minus_current, magnitude)
	current_fish["position"] = np.add(current_fish["position"], np.random.random() * step * direction)
	current_fish["fitness"] = objective(current_fish["position"])

def follow_state(current_fish, best_pos, step, objective):
	'''Executes the follow state behavior. Fish moves toward global best fish.

	Keyword arguments:
	pop (list) -- entire fish population, elements are dictionaries representing individual fish
	current_fish (dict) -- fish whose visual scope is calculated
	fish_in_vis (list) -- current number of main loop iteration
	try_number (int) -- number of times prey behavior is executed
	step (float) -- max step size fish can move
	objective (function) - objective function used to evaluate population fitness
	'''
	best_minus_current = np.subtract(best_pos, current_fish["position"])
	magnitude = np.absolute(np.subtract(best_pos, current_fish["position"]))
	direction = np.divide(best_minus_current, magnitude)
	current_fish["position"] = np.add(current_fish["position"], np.random.random() * step * direction)
	current_fish["fitness"] = objective(current_fish["position"])

def init_passed_population(population, pop_size, problem_size, visual, objective):
	"""Initializes population of fish based on an existing population.

	Keyword arguments:
	population (list) -- list of positions of all members of population
	pop_size (int) -- total number of bats in the population
	problem_size (int) -- dimension of problem (number of independent variables)
	visual (float) -- distance from fish to edge of visual (radius, if circular)
	objective (function) -- objective function used to evaluate population fitness

	Returns:
	list -- fish population, elements are dictionaries representing individual fish
	"""
	pop = np.array([{"position": population[i]} for i in range(pop_size)])
	for p in pop:
		p["fitness"] = objective(p["position"])
		# p["visual"] = calc_visual_scope(p, visual, problem_size)
	return pop

def search(objective, search_space, max_gens, pop,
		   visual=0.1, crowd=2, step=0.1, try_number=1):
	"""Performs artificial fish swarm algorithm search to find global minimum of passed objective function.

	Keyword arguments:
	objective (function) -- objective function to be minimized
	search_space (list) -- bounds for each dimension of the search space
	max_gens (int) -- maximum number of generations (iterations) the main loop will run for
	pop (list) -- list of positions of all initial members of population
	pop_count (int) -- total number of bats in the population
	visual (float) -- distance from fish to edge of visual (radius, if circular) (default = 0.1)
	crowd (int) -- number of fish that constitutes a crowd, used to determine behaviors

	Returns:
	list -- final bat positions
	"""
	pop_count = len(pop)
	crowd = int(round(crowd, 0))
	try_number = int(round(try_number, 0))
	problem_size = len(search_space)    # search space provides bounds for each dimension of problem,
										# length of this list provides number of dimensions
	fish = init_passed_population(pop, pop_count, problem_size, visual, objective)    # initialize passed population as fish
	best = min(fish, key=lambda x:x["fitness"])         # store intial best bat, based on lowest fitness
	for t in range(max_gens):
		for i in range(pop_count):      # foreach Fishi do
			calc_visual_scope(fish[i], visual, problem_size, search_space)         # Compute the visual scope
			fish_in_visual = which_fish_in_visual(fish, i, problem_size)
			num_in_visual = len(fish_in_visual)
			if (num_in_visual == 0):          # if Visual Scope is empty then
				move_state(fish[i], search_space, problem_size, objective)     # Y i   random(Xi)
			else:
				if (num_in_visual >= crowd):    # if Visual Scope is crowded then
					prey_state(fish, fish[i], fish_in_visual, try_number, step, objective) # Y i   search(Xi)
				else:
					center_position = calc_centroid(fish, fish_in_visual, num_in_visual, problem_size)        # Calculate Swarm centroid Ci 
					if (objective(center_position) < fish[i]["fitness"]):    # if Fitness of Ci better than Fitness of Fishi then
						swarm_state(fish[i], center_position, step, objective)     # Y i1   swarm(Xi)
					else:
						prey_state(fish, fish[i], fish_in_visual, try_number, step, objective)    # Y i1   search(Xi)
					if (fish[i]["fitness"] > best["fitness"]):   # if Fitness of Fishi worse than global best then
						follow_state(fish[i], best["position"], step, objective)     # Y i2   chase(Xi)
					else:
						prey_state(fish, fish[i], fish_in_visual, try_number, step, objective)    # Y i2   search(Xi)
					#y = min(y1, y2)          # Y i = argminfitness(Y i1 ; Y i2 )
			# check for best here
			if (fish[i]["fitness"] < best["fitness"]):
				fish[i]["position"] = best["position"]
				fish[i]["fitness"] = best["fitness"]
	final_positions = [fish[i]["position"] for i in range(pop_count)]
	if __name__ == "__main__":
		print("best =", best["position"], "fitness =", best["fitness"])    # un-comment to print out results
	return final_positions

	#######################
	# WILL IMPLEMENT LAST #
	#######################
	#     foreach Fishi do
	#         Xi   select(Xi; Y i)
	#         if No Advance In global best for T Iterations then
	#             Randomly choose a fish Xi Y i   leap(Xi)

def main():
	# initialize parameters
	pop_size = 25           # population size
	max_iterations = 100
	#tol = 0.00001           # stop tolerance (unused in this implementation)
	problem_size = 2        # dimensions of search space
	visual = 0.1
	crowd = 3
	step = 0.1
	try_number = 1

	search_space = [[-5.0, 5.0], [-5.0, 5.0]]

	initial_pop = [[np.random.uniform(search_space[j][0], search_space[j][1]) for j in range(problem_size)] for i in range(pop_size)]

	final_positions = search(objective_function, search_space, max_iterations, initial_pop, pop_size,
		visual, crowd, step, try_number)

if __name__== "__main__":
	main()

# def rand_visual(vis, problem_size):
#     rand_vis = [np.random.uniform(vis[i][0], vis[i][1]) for i in range(len(problem_size))]
#     return rand_vis

# def rand_move_in_visual(f, problem_size):
#     new_pos = np.zeros(problem_size)
#     for i in range(problem_size):
#         new_pos[i] = np.random.uniform(f["visual"][i][0], f["visual"][i][1])
#     return new_pos

# from paper:
# def move_state():
#     pass

# def prey_state(try_number, visual, step, objective):
#     for i in range(try_number):
#         # xj = xi + random(visual)
#         Yj = objective(xj)
#         # less than, to minimize
#         if (Yj < Yi):
#             xi_next = xi + random(step)*((xj-xi)/math.abs(xj-xi))
#         else:
#             xi_next = xi + random(step)
#     return objective(xi_next)

# def swarm_state(friend_number, visual, congestion):
#     friends = 0
#     xc = 0
#     for j in range(friend_number):
#         if (distance_ij < visual):
#             friends += 1
#             xc += xj
#     xc = xc/nf
#     # less than, to minimize
#     if (Yc/nf < congestion*Yi):
#         xi_next = xi + random(step)*((xc-xi)/math.abs(xc-xi))
#     else:
#         prey_state()
#     return objective(xi_next)

# def follow_state(friend_number, visual, congestion):
#     Ymin = Yj
#     for j in range(friend_number):
#         if (distance_ij < visual and Yj < Ymin):
#             Ymin = Yj       # fitness
#             xmin = xj       # position
#     friends = 0
#     for j in range(friend_number):
#         if (distance_minj < visual):
#             friends += 1
#     if (Ymin/friends < congestion*Yi):
#         xi_next = xi + random(step)*((xmin-xi)/math.abs(xmin-xi))
#     else:
#         prey_state()
#     return objective(xi_next)