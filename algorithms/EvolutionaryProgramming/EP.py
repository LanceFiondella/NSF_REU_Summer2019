# Evolutionary Programming algorithm in the Python Programming Language

# Adapted from:
# The Clever Algorithms Project: http://www.CleverAlgorithms.com
# (c) Copyright 2010 Jason Brownlee. Some Rights Reserved.
# This work is licensed under a Creative Commons Attribution-Noncommercial-Share Alike 2.5 Australia License.

import math
import random
from random import randint
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random


def fun(x, y):
    return x ** 2 + y ** 2


def objective_function(vector):                                                                                         # Calculate Z value (How 'high' from the min) This function contains what the minimums look like
    return ((vector[1]**2.0) + (vector[0]**2.0))



def random_vector(minmax):                                                                                              # Gives random numbers between the bounds set [-5, 5] for random vector and [0, 0.5] for strategy
    rand_vec = []
    for i in range(0, len(minmax)):
        rand_vec.append(minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random.random()))
    return rand_vec



def random_gaussian(mean=0.0, stdev=1.0):                                                                               # Standard normal distribution that is from [-3, 3]
  u1 = u2 = w = 0
  while True:
    u1 = 2 * random.random() - 1
    u2 = 2 * random.random() - 1
    w = u1 * u1 + u2 * u2
    if not (w >= 1):
        break
  w = math.sqrt((-2.0 * math.log(w)) / w)
  return (mean + (u2 * w) * stdev)




def mutate(candidate, search_space):
    child = {'vector':[], 'strategy':[]}
    for i in range(0, len(candidate['vector'])):
        s_old = candidate['strategy'][i]
        v = candidate['vector'][i] + s_old * random_gaussian()
        if v < search_space[i][0]: v = search_space[i][0]
        if v > search_space[i][1]: v = search_space[i][1]
        child['vector'].append(v)
        child['strategy'].append(s_old + random_gaussian() * abs(s_old) ** 0.5)
    return child



def tournament(candidate, population, bout_size):
  candidate['wins'] = 0
  for i in range(0, bout_size):
    other = population[randint(0,len(population)-1)]
    if candidate['fitness'] < other['fitness']: candidate['wins'] += 1



def init_population(minmax, pop_size):
    strategy = []
    for i in range(0, len(minmax)):
        strategy.append([0,  (minmax[i][1]-minmax[i][0]) * 0.05])
    pop = [dict() for x in range(0, pop_size)]
    for i in range(0, len(pop)):
        pop[i]['vector'] = random_vector(minmax)
        pop[i]['strategy'] = random_vector(strategy)
    for i in pop:
        i['fitness']= objective_function(i['vector'])
    return pop



def search (max_gens, search_space, pop_size, bout_size):
    population = init_population(search_space, pop_size)
    for i in population:
        i['fitness'] = objective_function(i['vector'])
    best = sorted(population, key=lambda i: i['fitness'])[0]

    for gen in range(0, max_gens):
        best = sorted(population, key=lambda i: i['fitness'])[0]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x = y = np.arange(-5, 5, 0.5)
        X, Y = np.meshgrid(x, y)
        zs = np.array([fun(x, y) for x, y in zip(np.ravel(X), np.ravel(Y))])
        Z = zs.reshape(X.shape)

        ax.plot_wireframe(X, Y, Z)

        ax.set_xlabel('X ')
        ax.set_ylabel('Y ')
        ax.set_zlabel('Z ')
        ax.view_init(70, 90)
        plt.plot([0], [0], [0], 'ro')
        for i in population:
            plt.plot([i['vector'][0]], [i['vector'][1]], [i['fitness']], 'go')
        plt.savefig(str(gen) + '.png')


        children = []
        for i in range(0,pop_size): children.append(mutate(population[i], search_space))
        for c in children: c['fitness'] = objective_function(c['vector'])
        children = sorted(children, key=lambda i: i['fitness'])
        if children[0]['fitness'] < best['fitness']: best = children[0]
        union = children + population
        for c in union: tournament(c, union, bout_size)
        union = sorted(union, key=lambda i: i['wins'])[::-1]
        population = union[0:pop_size]
        print (" > gen " + str(gen) + " fitness = " + str(best['fitness']) + " vector = " + str(best['vector'] ) + " strat = " + str(best['strategy'] ))


    return best



if __name__== "__main__":
    # problem configuration
    problem_size = 2
    search_space = []
    for i in range(0,problem_size): search_space.append([-5, +5])
    max_gens = 200
    pop_size = 100
    bout_size = 5
    # execute the algorithm
    best = search(max_gens, search_space, pop_size, bout_size)
    print(" done! Solution: f= " + str(best['fitness']) + " s = " + str(best['vector']))




