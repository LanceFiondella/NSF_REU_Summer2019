# Greedy Randomized Adaptive Search Procedure in the Ruby Programming Language

# The Clever Algorithms Project: http://www.CleverAlgorithms.com
# (c) Copyright 2010 Jason Brownlee. Some Rights Reserved.
# This work is licensed under a Creative Commons Attribution-Noncommercial-Share Alike 2.5 Australia License.

def euc_2d(c1, c2)
  Math.sqrt((c1[0] - c2[0])**2.0 + (c1[1] - c2[1])**2.0).round
end

def cost(perm, cities)
  distance =0
  perm.each_with_index do |c1, i|
    c2 = (i==perm.size-1) ? perm[0] : perm[i+1]
    distance += euc_2d(cities[c1], cities[c2])
  end
  return distance
end

def stochastic_two_opt(permutation)
  perm = Array.new(permutation)
  c1, c2 = rand(perm.size), rand(perm.size)
  exclude = [c1]
  exclude << ((c1==0) ? perm.size-1 : c1-1)
  exclude << ((c1==perm.size-1) ? 0 : c1+1)
  c2 = rand(perm.size) while exclude.include?(c2)
  c1, c2 = c2, c1 if c2 < c1
  perm[c1...c2] = perm[c1...c2].reverse
  return perm
end

def local_search(best, cities, max_no_improv)                                              #THE SEARCH CAN NOW COMMENCE WITH (BEST = HASH TABLE OF CANDIDATES AND THAT OPTIMAL COST), (CITIES = ORIGINAL DATA ARRAY), (MAX_NO_IMPROV = NUMBER OF TIMES CHECKED THAT THE COST DOES NOT IMPROVE)
  count = 0                                                                                #THIS MONITORS HOW MANY TIMES THE BEST IS STILL...THE BEST. AFTER 50 TIMES IN THIS CASE, OF THE BEST BEING THE BEST, WITH NO OTHER CANDIDATES TAKING IT'S PLACE, THE SEARCH WILL STOP
  begin
    candidate = {:vector=>stochastic_two_opt(best[:vector])}                               #FIRST, A RANDOM VECTOR IS CREATED WITH, THE ORIGINAL "GREEDY" VECTOR
    candidate[:cost] = cost(candidate[:vector], cities)                                    #THE NEW CONTENDER HAS IT'S COST CALCULATED AGAIN
    count = (candidate[:cost] < best[:cost]) ? 0 : count+1                                 #THEN A COMPARISION IS MADE TO SEE IF THE COST OF THE CANDIDATE IS LESS THAN THE CURRENT BEST (LOWER IS BETTER)
    best = candidate if candidate[:cost] < best[:cost]                                     #IF SO, THE CANDIDATE COST THEN TAKES THE PLACE OF THE BEST COST AND THE COUNT RESETS TO 0 TO START THE LOOPS OVER
  end until count >= max_no_improv                                                         #BASICALLY, THE BEST MUST REMAIN THE BEST FOR A SET TIME AND IN THIS CASE, THE EXAMPLE CHOOSES 50
  return best
end

def construct_randomized_greedy_solution(cities, alpha)                                   #THIS CREATES A INITIAL SAMPLE
  candidate = {}                                                                          #CREATES HASH THAT WILL HOLD VECTOR OF DIFFERENT "CITIES"
  candidate[:vector] = [rand(cities.size)]                                                #CHOSES RANDOM INDEX TO START THE VECTOR IN CANDIDATE HASH
  allCities = Array.new(cities.size) {|i| i}                                              #CREATES A NEW ARRAY WITH SIZE OF THE ORIGINAL CITIES ARRAY THAT CONTAINS ALL INDICES (0-51 IN THIS CASE)
  while candidate[:vector].size < cities.size                                             #THIS WHILE LOOP WITH CONTINUE AS THE CANDIDATE VECTOR FILLS UP
    candidates = allCities - candidate[:vector]                                           #TAKES THE ELEMENTS RANDOM ELEMENTS THAT ARE IN THE CANDIDATE VECTOR OUT OF ALL-CITIES AND PUTS THAT VALUE INTO A VECTOR CALLED CANDIDATES (THUS AS CANDIDATE VECTOR FILLS UP, CANDIDATES GETS SMALLER)
    costs = Array.new(candidates.size) do |i|                                             #DOES CALCULATIONS FOR CANDIDATES AND PUTS THEM INTO A COSTS ARRAY
      euc_2d(cities[candidate[:vector].last], cities[i])
    end
    rcl, max, min = [], costs.max, costs.min                                              #CREATES A RCL, MAX, AND MIN VARIABLES THAT WILL BE USED TO GAUGE IF THEY BELONG IN THE RCL
    costs.each_with_index do |c,i|
      rcl << candidates[i] if c <= (min + alpha*(max-min))                                #IF THE COST
    end
    candidate[:vector] << rcl[rand(rcl.size)]                                             #CHOSES RANDOM VALUE INSIDE FROM THE RCL AND APPENDS TO THE CANDIDATE VECTOR.
  end
  candidate[:cost] = cost(candidate[:vector], cities)                                     #WHAT YOU END UP WITH IS A CANDIDATE VECTOR THAT ALSO HAS A BASE COST THAT WILL BE LATER REFINED IN THE SEARCH PROCESS
  return candidate
end

def search(cities, max_iter, max_no_improv, alpha)
  best = nil
  max_iter.times do |iter|
    candidate = construct_randomized_greedy_solution(cities, alpha);                        #THIS BUILDS A SORT OF BASELINE VECTOR AND COST ASSOCIATED WHICH IS TO BE COMPARED WITH DURING THE SEARCH PHASE
    candidate = local_search(candidate, cities, max_no_improv)                              #HERE, THE PREVIOUS BASELINE SAMPLE IS COMPARED WITH MULTIPLE VARIATIONS OF SOLUTIONS UNTIL A BETTER OR SAME SOLUTIONS IS FOUND THAT HAS SAME OR EQUAL DISTANCE
    best = candidate if best.nil? or candidate[:cost] < best[:cost]                         #IF THE CURRENT BEST IS NULL(FIRST TIME) OR THE CANDIDATE HAS A LOWER COST THAN THE CURRENT BEST THEN THAT CANDIDATE BECOMES THE NEW BEST
    puts " > iteration #{(iter+1)}, best=#{best[:cost]} v=#{best[:vector].inspect}"         #PRINTS THE CURRENT ITERATION
  end
  return best
end

if __FILE__ == $0
  # problem configuration
  berlin52 = [[565,575],[25,185],[345,750],[945,685],[845,655],
              [880,660],[25,230],[525,1000],[580,1175],[650,1130],[1605,620],
              [1220,580],[1465,200],[1530,5],[845,680],[725,370],[145,665],
              [415,635],[510,875],[560,365],[300,465],[520,585],[480,415],
              [835,625],[975,580],[1215,245],[1320,315],[1250,400],[660,180],
              [410,250],[420,555],[575,665],[1150,1160],[700,580],[685,595],
              [685,610],[770,610],[795,645],[720,635],[760,650],[475,960],
              [95,260],[875,920],[700,500],[555,815],[830,485],[1170,65],
              [830,610],[605,625],[595,360],[1340,725],[1740,245]]                          #MAIN DATASET
  # algorithm configuration
  max_iter = 50                                                                          #NUMBER OF ITERATIONS
  max_no_improv =50                                                                        #NUMBER OF IMPROVEMENTS BEFORE A BEST IS SELECTED
  greediness_factor = 0.3                                                                  #HOW PICKY THE RESTRICTED CANDIDATE LIST IS WHICH IS A LIST OF OPTIMAL SOLUTIONS
  # execute the algorithm
  best = search(berlin52, max_iter, max_no_improv, greediness_factor)
  puts "Done. Best Solution: c=#{best[:cost]}, v=#{best[:vector].inspect}"
end
