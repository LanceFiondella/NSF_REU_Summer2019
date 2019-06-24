# Iterated Local Search algorithm in python

# adapted from Jason Brownlee's The Clever Algorithms Project: http://www.CleverAlgorithms.com
#   by Jacob Aubertine, 22 May 2019
# random library used for generating random numbers
import random
import math

# Euclidean (straight-line) distance between two points
def euc_2d(c1, c2):
    distance = round(math.sqrt((c1[0] - c2[0])**2.0 + (c1[1] - c2[1])**2.0))
    return distance

def cost(permutation, cities):
    distance = 0
    # iterates over elements of permutation; i is index, c1 is value
    for i, c1 in enumerate(permutation):
        # if (i==permutation.size-1)
        #   c2 = permutation[0]
        # else
        #   c2 = permutation[i+1]
    
        # c1 starts as first element of list, c2 is next element of list
        #   only time this isn't true is for last element: c1 = last city, c2 = first city
        c2 = permutation[0] if (i == len(permutation)-1) else permutation[i+1]
        distance += euc_2d(cities[c1], cities[c2])  # calculate distance from city 1 to city 2
    return distance

def random_permutation(cities):
    perm = list(range(len(cities))) # perm is list that holds value from 0 to 51 (52 total elements = number of cities)
    for i in range(len(perm)):
        #r = random.randint(0, len(perm)-i) + i
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
    #perm[c1:c2] = perm[c1:c2][::-1]
    return perm

def local_search(best, cities, max_no_improv):
    count = 0
    # continue local searches (2 opt) until max_no_improv iterations between improvement
    while (count < max_no_improv):
        candidate = {
            "vector": stochastic_two_opt(best["vector"])
        }   # create candidate dictionary, "vector" key maps to result of 2 opt on previous best vector
        candidate["cost"] = cost(candidate["vector"], cities)   # cost of candidate = total distance between current candidate path
        count = 0 if (candidate["cost"] < best["cost"]) else count+1    # if candidate cost is lower than previous best, improvement was made
        if (candidate["cost"] < best["cost"]):
            best = candidate    # candidate with lower cost becomes current best 
    return best

def double_bridge_move(perm):
    # (4-opt) divide the whole permutation into four approx equal parts (a, b, c, d) and reconnect them to find new solution
    # a = [0...pos1], b = [pos1...pos2], c = [pos2...pos3], d = [pos3..perm.size]
    # originally: a --> b --> c --> d
    # after:      a --> d --> c --> b

    pos1 = 1 + random.randint(0, len(perm) / 4)  # get start and end points of segments
    pos2 = pos1 + 1 + random.randint(0, len(perm) / 4)
    pos3 = pos2 + 1 + random.randint(0, len(perm) / 4)
    p1 = perm[0:pos1] + perm[pos3:]
    p2 = perm[pos2:pos3] + perm[pos1:pos2]
    return (p1 + p2)

def perturbation(cities, best):
    candidate = {}  # create new candidate dictionary
    candidate["vector"] = double_bridge_move(best["vector"])  # performs double bridge move using current best permutation
    candidate["cost"] = cost(candidate["vector"], cities) # calculates cost after perturbation
    return candidate

def search(cities, max_iterations, max_no_improv):
    best = {}
    best["vector"] = random_permutation(cities)  # vector is one key of dictionary, initialize to random permutation of cities
    best["cost"] = cost(best["vector"], cities) # cost is other key of dictionary
    best = local_search(best, cities, max_no_improv)  # initial local search provides new best
    for i in range(max_iterations):
        candidate = perturbation(cities, best)  # perform double bridge move as perturbation
        candidate = local_search(candidate, cities, max_no_improv)  # perform stochastic 2 opt for embedded local search
        if (candidate["cost"] < best["cost"]):
            best = candidate    # lowest cost (distance) is better, recorded
        print(" > iteration ", i+1, ", best = ", best["cost"], sep='')
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
    max_iterations = 100  # used in search method, does max_iterations local searches
    max_no_improv = 50  # used in local searches
    # execute the algorithm
    best = search(berlin52, max_iterations, max_no_improv)
    # print best solution cost, and vector of cities
    print("Done. Best solution: c =", best["cost"], ", v = [", end='')
    for x in range(len(best["vector"])-1):
        print(best["vector"][x], end=', ')
    print(best["vector"][-1], "]", sep='')

if __name__== "__main__":
  main()