# The following code executes the python version of the
# Random Search Function originally written in Ruby by Jason Brownlee
# newer version created 5/22 at 2pm
from random import randrange, uniform
import matplotlib.pyplot as plt
iteration_set=[]

def objective_function(vector):
    return sum([x**2 for x in vector])


def random_vector(mmx):
    return [uniform(mmx[x][0], mmx[x][1]) for x in range(len(mmx))]


def search(search_space, max_iter):
    best = None

    for iteration in range(max_iter):
        candidate = {}
        candidate["vector"] = random_vector(search_space)
        candidate["cost"]=objective_function(candidate["vector"])
        if(best == None) or (candidate["cost"] < best["cost"]):
            best = candidate
        iteration_set.append(best["cost"])
        print("iteration= {},best={}".format( iteration+1, best["cost"]))
    return best


# problem configuration
problem_size = 2
search_space = [[-5, 5] for x in range(problem_size)]
# algorithm configuration
max_iter = 100
# execute the algorithm
best = search(search_space, max_iter)
print("Done. Best solution: cost = {}, vector = {}".format( best["cost"], best["vector"]))

print(iteration_set)

plt.plot(iteration_set)
plt.show()