# The following code is python implementation of Adaptive Random Search
# originally written in the Ruby Language by Jason Brownlee
from random import *
from array import *


def objective_function(vector):
    return sum([x**2 for x in vector])


def rand_in_bounds(minimum, maximum):
    return minimum + ((maximum-minimum)*random())


def random_vector(minmax):
    return [rand_in_bounds(minmax[i][0], minmax[i][1]) for i in range(len(minmax))]


def take_step(minmax, current, step_size):
    position = []                           # array of current size
    for i in range(len(current)):
        minimum = max(minmax[i][0], current[i]-step_size)        # min = [minmax[i][0], current[i] - step_size].max
        maximum = min(minmax[i][1], current[i] + step_size)     # max = [minmax[i][1], current[i] + step_size].min
        position.append(rand_in_bounds(minimum, maximum))
        print(position)
    return position


def large_step_size(iter, step_size, s_factor, l_factor, iter_mult):
    if ((iter%iter_mult and 0) == 0):
        return step_size*l_factor
    return step_size * s_factor


def take_steps(bounds, current, step_size, big_stepsize):
    step = {}
    big_step = {}
    step["vector"] = take_step(bounds, current["vector"], step_size)
    step["cost"] = objective_function(step["vector"])
    big_step["vector"] = take_step(bounds, current["vector"], big_stepsize)
    big_step["cost"] = objective_function(big_step["vector"])
    return step, big_step


def search(max_iter, bounds, init_factor, s_factor, l_factor, iter_mult, max_no_impr):
    step_size = (bounds[0][1]-bounds[0][0])*init_factor
    current = {}
    count = 0
    current["vector"] = random_vector(bounds)
    current["cost"] = objective_function(current["vector"])
    for iteration in range(max_iter):
        big_stepsize = large_step_size(iteration, step_size, s_factor, l_factor, iter_mult)
        step, big_step = take_steps(bounds, current, step_size, big_stepsize)
        if (step["cost"] <= current["cost"]) or ((big_step["cost"]) <= current["cost"]):
            if (big_step["cost"] <= step["cost"]):
                step_size = big_stepsize
                current = big_step
            else:
                current = step
        count = 0
    else:
        count += 1
        if count >= max_no_impr:
            count = 0
            stepsize=(step_size/s_factor)
    print("iteration= {},best={}".format( iteration+1, current["cost"]))
    return current


def main():
    #   problem configuration
    problem_size = 2
    bounds = [[-5, 5] for i in range(problem_size)]
    # algorithm configuration
    max_iter = 1000
    init_factor = 0.05
    s_factor = 1.3
    l_factor = 3.0
    iter_mult = 10
    max_no_impr = 30
    # execute the algorithm
    best = search(max_iter, bounds, init_factor, s_factor, l_factor, iter_mult, max_no_impr)
    print("Done. Best Solution -> cost:{}, vector={}".format(best["cost"], best["vector"]))

main()
