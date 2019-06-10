# Greedy Randomized Adaptive Search Procedure in the Python Programing Language

# Adapted from:
# The Clever Algorithms Project: http://www.CleverAlgorithms.com
# (c) Copyright 2010 Jason Brownlee. Some Rights Reserved.
# This work is licensed under a Creative Commons Attribution-Noncommercial-Share Alike 2.5 Australia License.
import math
import time
from random import randint
from matplotlib import pyplot as plt
from matplotlib import animation


# Global Arrays
iteration_set = []
run_time = []
changed_factor = []
bests = []
iters = []


# Distance Formula
def euc_2d(c1,c2):
    return round(math.sqrt((c1[0] - c2[0])**2.0 + (c1[1] - c2[1])**2.0))                                                # Formula to Calculate Distances Between Points


# Calculate Total Distance of Vector
def cost(perm, cities):
    distance = 0                                                                                                        # Initialize distance as zero
    for i in range(0, len(perm)):                                                                                       # Continues for lenght of vector
        if i == len(perm)-1:                                                                                            # If the vector has reached the end, the second endpoint is the first element
            c2 = perm[0]                                                                                                # Else the second point is the next element in perm
        else:
            c2 = perm[i+1]
        distance += euc_2d(cities[perm[i]], cities[c2])                                                                 # Keep adding to distance to find total distance
    return distance


# Calculate permutation of existing vector.
def stochastic_two_opt(permutation):
    perm = permutation.copy()                                                                                           # Makes copy of permutation
    c1, c2 = randint(0,len(perm)), randint(0,len(perm))                                                                 # Chooses two endpoints from the permutation(indicies)
    exclude = [c1]                                                                                                      # Creates exclude list that includes adjacent points before and after c1
    exclude.append(len(perm)-1) if c1 == 0 else exclude.append(c1-1)
    exclude.append(0) if c1 == len(perm)-1 else exclude.append(c1+1)
    while exclude.__contains__(c2) : c2 = randint(0,len(perm))                                                          # Makes sure c2 is not within the range of c1-1, c1, or c+1
    if c2 < c1 : c1, c2 = c2, c1                                                                                        # Makes sure c1 is less than c2
    perm[c1:c2] = perm[c1:c2][::-1]                                                                                     # Takes the section between c1 and c2, reverses it and inserts it back into the vector
    return perm


# Refining search for greedy solution
def local_search(best, cities, max_no_improv):
    count = 0                                                                                                           # Initializes count to 0
    while count <= max_no_improv:                                                                                       # Continues while the count is less than or equal to max_no_improv
        candidate = {'vector': stochastic_two_opt(best['vector'])}                                                      # Creates new random vector using Stochastic Two Opt method
        candidate['cost'] = cost(candidate['vector'], cities)                                                           # Calculates total cost of vector
        if candidate['cost'] < best['cost']: count = 0                                                                  # If the new cost is better than the current best, the counter starts over, else the counter increments
        else: count = count+1
        if candidate['cost'] < best['cost']: best = candidate                                                           # If the counter starts over, meaning there is a new best, the candidate becomes new best
    return best


# Creates first greedy solution
def construct_randomized_greedy_solution(cities, alpha):
    candidate = {'vector': [randint(0,len(cities)-1)]}                                                                  # Starts the greedy solution vector with a random city
    allCities = range(0,len(cities))                                                                                    # Creates a list of all cities (All indices)
    while len(candidate['vector']) < len(cities):                                                                       # Continues while the greedy solution vector is not done (All cities should be represented)
        candidates = list(set(allCities) - set(candidate['vector']))                                                    # Cadidates only holds the remaining cities to be put into the solutions vector
        costs=[]
        for i in range(0,len(candidates)):
            costs.append(euc_2d(cities[candidate['vector'][-1]], cities[i]))                                            # Append all distances from the last element in the solution vector to all remaining points that need to be put into the solution
        rcl, maximum, minimum = [], max(costs), min(costs)                                                              # Starts the selection progress for which points qualify in the greedy selection process
        for i in range(0, len(costs)):                                                                                  # Sorts through all the distances, and if they point qualifies, it is put into a restricted candidate list
            if costs[i] <= (minimum + alpha *(maximum-minimum)): rcl.append(candidates[i])
        candidate['vector'].append(rcl[randint(0,len(rcl)-1)])                                                          # Append a random value from that RCL to be the next city in the greedy solution vector
    candidate['cost'] = cost(candidate['vector'], cities)                                                               # Finally calculate the total cost (distance) of that greedy solution vector
    return candidate


# Main search routine
def search(cities, max_iter, max_no_improv, alpha):
    best = None                                                                                                         # The current best starts off as None
    for i in range(0, max_iter):                                                                                        # Iterate for a set number of iterations set by max_iter
        candidate = construct_randomized_greedy_solution(cities, alpha)                                                 # Create a greedy solution
        candidate = local_search(candidate, cities, max_no_improv)                                                      # Refine search of the solution
        if best == None or candidate['cost'] < best['cost']: best = candidate                                           # If there is a new best, replace the old best
        iteration_set.append(best['vector'])                                                                            # Vectors for graphing
        iters.append(i+1)
        #bests.append(best['cost'])
        print("> iteration " + str(i+1) + " best = " + str(best['cost']) + " vector = " + str(best['vector']))
    return best


# Main function
def main(max_iter, max_no_improv, greediness_factor):
    berlin52 =[[565, 575], [25,185], [345,750], [945,685], [845, 655],
                [880, 660], [25, 230], [525, 1000], [580, 1175], [650, 1130], [1605, 620],
                [1220, 580], [1465, 200], [1530, 5], [845, 680], [725, 370], [145, 665],
                [415, 635], [510, 875], [560, 365], [300, 465], [520, 585], [480, 415],
                [835, 625], [975, 580], [1215, 245], [1320, 315], [1250, 400], [660, 180],
                [410, 250], [420, 555], [575, 665], [1150, 1160], [700, 580], [685, 595],
                [685, 610], [770, 610], [795, 645], [720, 635], [760, 650], [475, 960],
                [95, 260], [875, 920], [700, 500], [555, 815], [830, 485], [1170, 65],
                [830, 610], [605, 625], [595, 360], [1340, 725], [1740, 245]]

    start = time.time()
    best = search(berlin52, max_iter, max_no_improv, greediness_factor)
    elapsed_time = (time.time() - start)
    print("Done. Best Solution " + str(best['cost']) + " vector = " + str(best['vector']))
    print("Time = " + str(elapsed_time))
    bests.append(best['cost'])
    return elapsed_time


if __name__== "__main__":

    x, y, z = 50, 50, 0.01
    #main(x,y,z)

    for i in range(0, 100):
        run_time.append(main(x,y,z))
        changed_factor.append(z)
        z+=0.01
    '''
    bests = bests[:-1]
    plt.plot(iters, bests)
    plt.ylabel('Cost')
    plt.xlabel('Number of Iterations')
    plt.title('Cost vs Number of Iterations')
    plt.show()
    '''
    plt.plot(changed_factor, run_time)
    plt.ylabel('Run Time')
    plt.xlabel('Value of Greediness')
    plt.title('Run Time to Decreasing Greediness')
    plt.show()

    plt.plot(changed_factor, bests)
    plt.ylabel('Best Cost')
    plt.xlabel('Value of Greediness')
    plt.title('Best Cost to Decreasing Greediness')
    plt.show()
dataset =[[565,575],[25,185],[345,750],[945,685],[845,655],
            [880,660],[25,230],[525,1000],[580,1175],[650,1130],[1605,620],
            [1220,580],[1465,200],[1530,5],[845,680],[725,370],[145,665],
            [415,635],[510,875],[560,365],[300,465],[520,585],[480,415],
            [835,625],[975,580],[1215,245],[1320,315],[1250,400],[660,180],
            [410,250],[420,555],[575,665],[1150,1160],[700,580],[685,595],
            [685,610],[770,610],[795,645],[720,635],[760,650],[475,960],
            [95,260],[875,920],[700,500],[555,815],[830,485],[1170,65],
            [830,610],[605,625],[595,360],[1340,725],[1740,245]]

iterations = 100

# Animation Methods by Josh Steakelum
def animate(i):							                                                    # function is ran every frame with the frame number as the argument
    if(i >= len(iteration_set)):		                                                    # do nothing if the number is over the number of frame data lists
        return lines
    x, y = [dataset[x][0] for x in iteration_set[i]], [dataset[x][1] for x in iteration_set[i]]
    for line in lines:					                                                    # look up the new line positioning in the saved frame data
        line.set_data(x,y) 				                                                    # assign the values to the lines in the graph
    return lines

def findaxes(dset): 					                                                    # supply this with the set of actual data values
    xs = [dset[x][0] for x in range(len(dset))]                                             # all it does is calculate the min/max and return some fitting rectangle around it
    ys = [dset[x][1] for x in range(len(dset))]
    xs.sort()
    ys.sort()
    return [	int(1.1*xs[0] - 0.1*xs[-1]), int(1.1*xs[-1] - 0.1*xs[0]),
			int(1.1*ys[0] - 0.1*ys[-1]), int(1.1*ys[-1] - 0.1*ys[0])]

fig, ax = plt.subplots()			                                                        # create plot objects
ax = plt.axis(findaxes(dataset)) 	                                                        # set the viewing area to include all data points

x, y = [dataset[x][0] for x in iteration_set[0]], [dataset[x][1] for x in iteration_set[0]] # find initial dot/line positions
dots = plt.plot(x,y,'o')			                                                        # plot them - dots don't move, so this one function call is sufficient
lines = plt.plot(x,y)

anim = animation.FuncAnimation(fig, animate, frames=range(int(iterations*1.5)), interval = 10, repeat=True)

plt.ylabel('Y')
plt.xlabel('X')
plt.title('Max Iter = 10000, Max No Improv = 1000, Greediness = 0.3')
anim.save('animation.gif', writer='imagemagick', fps=60)
plt.show()
