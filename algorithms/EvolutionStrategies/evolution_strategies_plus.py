# Evolution Strategies algorithm in python

# adapted from Jason Brownlee's The Clever Algorithms Project: http://www.CleverAlgorithms.com
#   by Jacob Aubertine, 30 May 2019

# random library used for generating random numbers, matplotlib for animation
import random
import math
from operator import itemgetter     # used as part of sorting lists of dictionaries by value
#from matplotlib import pyplot as plt
#from matplotlib import animation

def objective_function(vector):
    # objective function --> f = x1^2 + x2^2

    # probably don't need to be this elaborate with it
    # this was trying to do it like the Ruby code, but that uses inject() method
    # no direct analog in python, and the code below looks inefficient
    ''' 
    function_vals = []
    for i in range(len(vector)):
        dimension_sum = 0.0
        for j in range(len(vector)):
            dimension_sum += (vector[i][j] ** 2.0)
        function_vals.append(dimension_sum)
    return function_vals
    '''
    # objective function hard coded, but I would think in most cases it needs to be
    return (vector[0] ** 2.0 + vector[1] ** 2.0)

def random_vector(minmax):
    # returns random vector (values) to be evaluated by objective function
    rand_vect = [(minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random.random())) for i in range(len(minmax))]
    return rand_vect

def random_gaussian(mean=0.0, stdev=1.0):
    u1, u2, w = 0, 0, 0
    while True:
        u1 = 2 * random.random() - 1
        u2 = 2 * random.random() - 1
        w = u1 * u1 + u2 * u2
        if not (w >= 1):
            break
    w = math.sqrt((-2.0 * math.log(w)) / w)
    return mean + (u2 * w) * stdev  # 99.72% of values within [-3, 3]


def mutate_problem(vector, stdevs, search_space):
    # vector is parent vector, stdevs is parent strategy, search_space is minmax
    child = []
    for i, v in enumerate(vector):
        child.append(v + stdevs[i] * random_gaussian())    # create child vector using parent vector + strategy * random std normal dist value
        # ensure new values within search space
        if (child[i] < search_space[i][0]):         
            child[i] = search_space[i][0]               # set value to lower bound if it was below it (-5)
        if (child[i] > search_space[i][1]):
            child[i] = search_space[i][1]               # set value to upper bound if it was above it (+5)
    return child

def mutate_strategy(stdevs):
    # mutate parent strategy (stdevs) to create child strategy
    # don't quite understand WHY this is the math yet, will look into sources beyond Clever Algorithms
    tau = math.sqrt(2.0 * float(len(stdevs))) ** -1.0
    tau_p = math.sqrt(2.0 * math.sqrt(float(len(stdevs)))) ** -1.0
    # create new child strategy list based on mutation of parent strategy
    child = [(stdevs[i] * math.exp(tau_p * random_gaussian() + tau * random_gaussian())) for i in range(len(stdevs))]
    return child

def mutate(par, minmax):
    # population[i] is par (parent), minmax is search_space x:(-5, 5), y:(-5, 5)
    child = {
        "vector": mutate_problem(par["vector"], par["strategy"], minmax),   # create child vector by mutating parent vector
        "strategy": mutate_strategy(par["strategy"])                        # create child strategy by mutating parent strategy
    }
    return child

def init_population(minmax, pop_size):
    # generate initial population
    strategy = [[0, (minmax[i][1] - minmax[i][0]) * 0.05] for i in range(len(minmax))]  
        # search spaced passed as minmax, (5 - (-5) * 0.05) = 10 * 0.05 = 0.5
        # strategy becomes list of lists = [[0, 0.5], [0, 0.5]]
        # referred to as "standard deviations"
    pop = [{} for i in range(pop_size)]    # pop initialized as empty list, will contain dictionary for each of (pop_size) elements
    for i in range(pop_size):
        pop[i]["vector"] = random_vector(minmax)                    # vector based on search space, [x,y] values
        pop[i]["strategy"] = random_vector(strategy)                # vector based on strategy
        pop[i]["fitness"] = objective_function(pop[i]["vector"])    # calculates fitness value for each member of population
                                                                    # firness is value of f when passed vector [x, y]
    return pop

def search(max_gens, search_space, pop_size, num_children):
    # main search loop
    population = init_population(search_space, pop_size)
    
    # need to sort values of all dictionaries contained in list
    best = sorted(population, key=itemgetter("fitness"))[0]     # want to minimize function
                                                                # sort members by fitness value lowest to highest
                                                                # best member of population is one with lowest fitness value
                                                                # apparently, itemgetter() quicker than lambda function

    for gen in range(max_gens):       # main search loop
        children = []    # children initialized as empty list, will contain dictionary for each of (num_children) elements
        for i in range(num_children):   
            children.append(mutate(population[i], search_space))   # vector and strategy for children mutated from their parents' vector and strategy
            children[i]["fitness"] = objective_function(children[i]["vector"])  # evaluate child fitness based on objective function
        union = sorted(children + population, key=itemgetter("fitness"))        # create union of children and parents (concatenate arrays)
        #union.sort(population, key=itemgetter("fitness"))   # sort children and parents by firness lowest to highest (best to worst)
        if (union[0]["fitness"] < best["fitness"]):         
            best = union[0]                                 # if new member with lowest fitness is better than previous best, it becomes new best
        population = union[:pop_size]                       # new population (next parents) become best 30 of 50 members
        print(" > gen", gen, "fitness =", best["fitness"])
    return best

def main():
    # problem configuration
    problem_size = 2        # need to test if this does more than make objective function x^2 + y^2
    search_space = [[-5, 5] for i in range(problem_size)]   # creates two lists = [-5, 5], x and y bounds
    # algorithm configuration
    max_gens = 100          # iterations for main loop, used as stop condition
    pop_size = 30           # mu, number of parents
    num_children = 20       # lambda, number of children
    # execute the algorithm
    best = search(max_gens, search_space, pop_size, num_children)     # main search loop
    print("Done! Solution: f =", best["fitness"], ", s = [", end='')
    for x in range(len(best["vector"])-1):
        print(best["vector"][x], end=', ')
    print(best["vector"][-1], "]", sep='')

if __name__== "__main__":
  main()