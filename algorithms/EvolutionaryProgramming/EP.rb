# Evolutionary Programming algorithm in the Ruby Programming Language

# The Clever Algorithms Project: http://www.CleverAlgorithms.com
# (c) Copyright 2010 Jason Brownlee. Some Rights Reserved.
# This work is licensed under a Creative Commons Attribution-Noncommercial-Share Alike 2.5 Australia License.

def objective_function(vector)
  return vector.inject(0.0) {|sum, x| sum +  (x ** 2.0)}                                                                # Calculates fitness value. sums up vector[0]^2 + vector[1]^2 and that's the fitness
end

def random_vector(minmax)
  return Array.new(minmax.size) do |i|
    minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * rand())
  end
end

def random_gaussian(mean=0.0, stdev=1.0)
  u1 = u2 = w = 0
  begin
    u1 = 2 * rand() - 1                                                                                                 # rand() generates float from 0 to 1.0
    u2 = 2 * rand() - 1
    w = u1 * u1 + u2 * u2                                                                                               # u1^2 + u2^2 = w
  end while w >= 1                                                                                                      # Make sure w is less than or equal to 1
  w = Math.sqrt((-2.0 * Math.log(w)) / w)                                                                               # Do calculations to w
  return mean + (u2 * w) * stdev                                                                                        # Return w with mean and stdev parameters
end

def mutate(candidate, search_space)
  child = {:vector=>[], :strategy=>[]}                                                                                  # Creates empty hash called child with 'vector' key and 'strategy' key that both have empty arrays
  candidate[:vector].each_with_index do |v_old, i|                                                                      # Iterate over candidate[:vector] which has 2 elements, v_old = value, i = index
    s_old = candidate[:strategy][i]                                                                                     # s_old = strategy old
    v = v_old + s_old * random_gaussian()                                                                               # v = strategy old + vector old * random_gaussian()
    v = search_space[i][0] if v < search_space[i][0]                                                                    # if v is less than minimum search space, v equals the minimum search space. (Contains v within)
    v = search_space[i][1] if v > search_space[i][1]                                                                    # if v is greater than maximum search space, v equals the maximum search space.
    child[:vector] << v                                                                                                 # Appends v to child hash under lab 'vector'
    child[:strategy] << s_old + random_gaussian() * s_old.abs**0.5                                                      # Appends a calculation including random_gaussian() to the strategy
  end
  return child                                                                                                          # Returns child hash
end

def tournament(candidate, population, bout_size)
  candidate[:wins] = 0                                                                                                  # Candidate is that single hash, put a new key value pair called 'wins' and initialize at 0
  bout_size.times do |i|                                                                                                # Do this loop for bout_size, in this case 5
    other = population[rand(population.size)]                                                                           # From the whole population, now with both children and parents so 200 elements since the original size was 100, pick a random index that points to a single 3 element hash and name that hash as other
    candidate[:wins] += 1 if candidate[:fitness] < other[:fitness]                                                      # if the current single candidate hash being tested has a lower fitness than the random fitness, mark it as a win in the candidate 'wins' key.
  end                                                                                                                   # Thus, a candidate can have a max of 5 wins and min of 0
end

#Why are all pop the same?
def init_population(minmax, pop_size)
  strategy = Array.new(minmax.size) do |i|                                                                              # Make a new array called strategy that holds [0, (5 - (-5) * 0.5)] problem_size times. [[], []] Nested array
    [0,  (minmax[i][1]-minmax[i][0]) * 0.05]
  end
  pop = Array.new(pop_size, {})                                                                                         # Array of size pop, 100 in this case of empty hashes
  pop.each_index do |i|                                                                                                 # For each element in pop which is a [{}, {}]
    pop[i][:vector] = random_vector(minmax)                                                                             # Run random_vector on the original search space and put it in the pop[index][:vector]
    pop[i][:strategy] = random_vector(strategy)                                                                         # Run random_vector on the strategy created earlier in this method and append it to pop[index][:strategy]
  end
  pop.each{|c| c[:fitness] = objective_function(c[:vector])}                                                            # Calculates a fitness for each vector using objective_function
  return pop
end

def search(max_gens, search_space, pop_size, bout_size)
  population = init_population(search_space, pop_size)
  population.each{|c| c[:fitness] = objective_function(c[:vector])}                                                     # Calculates fitness of each vector....again?
  best = population.sort{|x,y| x[:fitness] <=> y[:fitness]}.first                                                       # Sorts by fitness and the best is the first element: It is a hash with 'vector' 'strategy' 'fitness'
  max_gens.times do |gen|
    children = Array.new(pop_size) {|i| mutate(population[i], search_space)}                                            # Creates mutations of current pop that fall within search space and witch each hash in the population
    children.each{|c| c[:fitness] = objective_function(c[:vector])}                                                     # Calculates fitness for all child nodes
    children.sort!{|x,y| x[:fitness] <=> y[:fitness]}                                                                   # Sorts all the hashes by fitness number in ascending order from 0 = lowest, pop_size = highest
    best = children.first if children.first[:fitness] < best[:fitness]                                                  # best is now equal to the first (lowest fitness) child node if it is lower than the current best fitness. Keep in mind that best is a single 3 element hash with 'vector' = [], 'strategy' = [] and 'fitness' = float
    union = children+population                                                                                         # Combines both population and current children
    union.each{|c| tournament(c, union, bout_size)}                                                                     # Passes single hash, whole array of hashes, and bout_size which is 5 here, this will happen 200 times here and returns that hash set with 'wins' key value pair which tells how many times a random fitness in union is less than the one in that set
    union.sort!{|x,y| y[:wins] <=> x[:wins]}                                                                            # Sort by number of wins from 0 = most, union.size = least
    population = union.first(pop_size)                                                                                  # The first pop_size of the union become the new population
    puts " > gen #{gen}, fitness=#{best[:fitness]}"                                                                     # Prints out once each gen , fitness should get less and less
  end
  return best
end

if __FILE__ == $0
  # problem configuration
  problem_size = 2                                                                                                      #
  search_space = Array.new(problem_size) {|i| [-5, +5]}                                                                 # Makes array of 2 element arrays with value [-5, 5]
  # algorithm configuration
  max_gens = 200
  pop_size = 100
  bout_size = 5
  # execute the algorithm
  best = search(max_gens, search_space, pop_size, bout_size)
  puts "done! Solution: f=#{best[:fitness]}, s=#{best[:vector].inspect}"
end
