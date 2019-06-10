# to run this code snippet - first put it at the end of your algorithm code:
# need some list of 2-d coordinate lists called dataset:
# 	e.g. dataset = [[1,1], [100,100], ... [x,y]]; this contains the coordinates for the "cities" in the TSP problem
# some list of lists called iteration_set, containing the sequence of cities to travel to:
# 	e.g. [2, 1, 4, 100, ..., n]

# this code should handle the rest, feel free to change the arguments of FuncAnimationm
# let me know if there's any issues


def animate(i):							# function is ran every frame with the frame number as the argument
	if(i >= len(iteration_set)):		# do nothing if the number is over the number of frame data lists
		return lines
	x, y = [dataset[x][0] for x in iteration_set[i]], [dataset[x][1] for x in iteration_set[i]]
	for line in lines:					# look up the new line positioning in the saved frame data
		line.set_data(x,y) 				# assign the values to the lines in the graph
	return lines

def findaxes(dset): 					# supply this with the set of actual data values
	xs = [dset[x][0] for x in range(len(dset))] # all it does is calculate the min/max and return some fitting rectangle around it
	ys = [dset[x][1] for x in range(len(dset))]
	xs.sort()
	ys.sort()
	return [	int(1.1*xs[0] - 0.1*xs[-1]), int(1.1*xs[-1] - 0.1*xs[0]),
			int(1.1*ys[0] - 0.1*ys[-1]), int(1.1*ys[-1] - 0.1*ys[0])]


fig, ax = plt.subplots()			# create plot objects
ax = plt.axis(findaxes(dataset)) 	# set the viewing area to include all data points

x, y = [dataset[x][0] for x in iteration_set[0]], [dataset[x][1] for x in iteration_set[0]] # find initial dot/line positions
dots = plt.plot(x,y,'o')			# plot them - dots don't move, so this one function call is sufficient
lines = plt.plot(x,y)


# create an animation with 1.5*the amount of frames to give a delay at the end, 10ms between frames (can be changed), repeating
anim = animation.FuncAnimation(fig, animate, frames=range(int(max_iterations*1.5)), interval = 10, repeat=True)
#anim.save('traveling_salesman.gif', fps=30, writer='imagemagick')	# save the animation as a gif, however this doesn't work on windows (need imagemagick)

plt.show()
