# Guided Local Search algorithm in python

# adapted from Jason Brownlee's The Clever Algorithms Project: http://www.CleverAlgorithms.com
#   by Jacob Aubertine, 23 May 2019

# random library used for generating random numbers, math for sqrt()
import random
import math

# Euclidean (straight-line) distance between two points
def euc_2d(c1, c2):
    distance = round(math.sqrt((c1[0] - c2[0])**2.0 + (c1[1] - c2[1])**2.0))
    return distance

def random_permutation(cities):
    perm = list(range(len(cities))) # perm is list that holds value from 0 to 51 (52 total elements = number of cities)
    for i in range(len(perm)):
        r = random.randint(0, len(perm)-1)  # selects random city
        perm[r], perm[i] = perm[i], perm[r] # swaps random city with current city iteration
    return perm # return list of random city numbers

def stochastic_two_opt(permutation):
    # take a route that crosses over itself and reorder it so that it does not
    # - A   B -             - A - B -
    #     X         ==>     
    # - C   D -             - C - D -
    
    perm = list(permutation)    # perm is copy of permutation
    c1, c2 = random.randint(0, len(perm)), random.randint(0, len(perm))   # random cities
    exclude = [c1]
    
    #if (c1 == 0):
    #    exclude.append(len(perm)-1)
    #else:
    #    exclude.append(c1-1)
    exclude.append(len(perm)-1 if (c1 == 0) else c1-1)  # appends previous city to list (final if c1 is first city in list)
    #exclude.append(len(perm)-1) if (c1 == 0) else exclude.append(c1-1)  
    
    #if (c1 == len(perm)-1):
    #    exclude.append(0)
    #else:
    #    exclude.append(c1+1)
    exclude.append(0 if (c1 == len(perm)-1) else c1+1)  # appends next city to list (first if c1 is last city in list)
    #exclude.append(0) if (c1 == len(perm)-1) else exclude.append(c1+1)  
    
    while (c2 in exclude):
        c2 = random.randint(0, len(perm))   # randomly choose a new city 2 if it one of the three excluded cities
    if (c2 < c1):
        c1, c2 = c2, c1 # ensures c1 < c2
    perm[c1:c2] = reversed(perm[c1:c2]) # reverses elements between city 1 and city 2
    return perm

def augmented_cost(permutation, penalties, cities, l):
    # parameter name lambda replaced with l, don't want lambda function
    distance, augmented = 0, 0
    # iterates over elements of permutation; i is index, c1 is value
    for i, c1 in enumerate(permutation):
        # if (i==permutation.size-1)
        #   c2 = permutation[0]
        # else
        #   c2 = permutation[i+1]
    
        # c1 starts as first element of list, c2 is next element of list
        #   only time this isn't true is for last element: c1 = last city, c2 = first city
        c2 = permutation[0] if (i == len(permutation)-1) else permutation[i+1]
        if (c2 < c1):
            c1, c2 = c2, c1
        d = euc_2d(cities[c1], cities[c2])
        distance += d  # calculate distance from city 1 to city 2
        augmented += d + (l * (penalties[c1][c2]))
    return [distance, augmented]

def cost(cand, penalties, cities, l):
    cost, acost = augmented_cost(cand["vector"], penalties, cities, l)
    # python dictionaries are mutable and original objects passed to function
    cand["cost"], cand["aug_cost"] = cost, acost    # edit cand dictionary that's passed to function

def local_search(current, cities, penalties, max_no_improv, l):
    cost(current, penalties, cities, l) # cost and aug_cost keys of current dict changed in cost function
    count = 0
    # perform 2-opt until no improvements for max_no_improv number of iterations
    while (count < max_no_improv):
        candidate = {"vector": stochastic_two_opt(current["vector"])}   # new candidate dictionary, vector is result of 2_opt on current vector
        cost(candidate, penalties, cities, l)   # calculate cost and augmented cost of new candidate
        count = 0 if (candidate["aug_cost"] < current["aug_cost"]) else count+1    # if candidate cost is lower than previous best, improvement was made
        if (candidate["aug_cost"] < current["aug_cost"]):
            current = candidate # lower cost candidate becomes current best 
    return current

def calculate_feature_utilities(penal, cities, permutation):
    utilities = [0 for i in range(len(permutation))]    # list of size len(permutation), elements initialized to 0
    for i, c1 in enumerate(permutation):
        c2 = permutation[0] if (i == len(permutation)-1) else permutation[i+1]  # get city in next index (city 2 = first element if c1 = last element)
        if (c2 < c1):
            c1, c2 = c2, c1
        utilities[i] = euc_2d(cities[c1], cities[c2]) / (1.0 + penal[c1][c2])   # utility formula: utility = cost div (1 + current penalty)
    return utilities

def update_penalties(penalties, cities, permutation, utilities):
    max_util = max(utilities)   # Penalties are only updated for those features in a locally optimal solution that maximize utility
    for i, c1 in enumerate(permutation):
        c2 = permutation[0] if (i == len(permutation)-1) else permutation[i+1]  # get city in next index (city 2 = first element if c1 = last element)
        if (c2 < c1):
            c1, c2 = c2, c1
        if (utilities[i] == max_util):
            penalties[c1][c2] += 1  # penalties updated by adding 1 to the penalty for the future (a counter)
    return penalties

def search(max_iterations, cities, max_no_improv, l):
    '''
    Can assign best to None (like how best = nil initially in Ruby code), but
    cannot perform check (current["cost"] < best["cost"]) because 'best' is
    unsubscriptable. Attempted solution is to assign best using cost function
    before for loop.
    '''

    current = {"vector": random_permutation(cities)}    # create current dict, vector maps to random permutation of cities
    best = current
    penalties = [[0 for i in range(len(cities))] for j in range(len(cities))]   # creates 2D list (len(cities) X len(cities)) initialized to 0
    cost(best, penalties, cities, l)    # get initial values for best, so it can be compared to
    for i in range(max_iterations):
        current = local_search(current, cities, penalties, max_no_improv, l)
        utilities = calculate_feature_utilities(penalties, cities, current["vector"])
        update_penalties(penalties, cities, current["vector"], utilities)
        if (current["cost"] < best["cost"]):
            best = current  # new best permutation if cost is lower
        print(" > iteration ", i+1, ", best = ", best["cost"], "aug = ", best["aug_cost"], sep='')
    return best
    

def main():
    random.seed(a=None) # set seed to use system time
    # problem configuration
    berlin52 = [[565,575],[25,185],[345,750],[945,685],[845,655],
        [880,660],[25,230],[525,1000],[580,1175],[650,1130],[1605,620],
        [1220,580],[1465,200],[1530,5],[845,680],[725,370],[145,665],
        [415,635],[510,875],[560,365],[300,465],[520,585],[480,415],
        [835,625],[975,580],[1215,245],[1320,315],[1250,400],[660,180],
        [410,250],[420,555],[575,665],[1150,1160],[700,580],[685,595],
        [685,610],[770,610],[795,645],[720,635],[760,650],[475,960],
        [95,260],[875,920],[700,500],[555,815],[830,485],[1170,65],
        [830,610],[605,625],[595,360],[1340,725],[1740,245]]
    # algorithm configuration
    max_iterations = 150  # used in search method, does max_iterations local searches
    max_no_improv = 20  # used in local searches
    alpha = 0.3 # α ∈ (0, 1] (around 0.3 for TSP and 2-opt)
    local_search_optima = 12000.0 # cost of a local optimum found by a local search (why is it 12000?)
    l = alpha * (local_search_optima/float(len(berlin52))) # λ = (α · cost(optima))/N
    # execute the algorithm
    best = search(max_iterations, berlin52, max_no_improv, l)
    # print best solution cost, and vector of cities
    print("Done. Best solution: c =", best["cost"], ", v = [", end='')
    for x in range(len(best["vector"])-1):
        print(best["vector"][x], end=', ')
    print(best["vector"][-1], "]", sep='')

if __name__== "__main__":
    main()