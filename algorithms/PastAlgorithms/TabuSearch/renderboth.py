import TabuSearch_Render
import ReactiveTabuSearch_Render

import matplotlib.pyplot as plt


fig, ax = plt.subplots()
ax.plot(TabuSearch_Render.costs)
ax.plot(ReactiveTabuSearch_Render.costs)
ax.legend(["Tabu Search", "Reactive Tabu Search"])
plt.xlabel("Iteration")
plt.ylabel("Cost")
plt.title("Travel Cost vs Iteration (lower better)")

fig2, ax2 = plt.subplots()
ax2.plot(TabuSearch_Render.time_between)
ax2.plot(ReactiveTabuSearch_Render.time_between)
ax2.legend(["Tabu Search", "Reactive Tabu Search"])
plt.xlabel("Iteration")
plt.ylabel("Time")
plt.title("Iteration Number vs Iteration Time (lower better)")

plt.show()