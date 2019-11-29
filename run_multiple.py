#!/usr/bin/python3

import subprocess as sp

for i in [2,5,10,15,20,25]:
	sp.run(["python3", "evaluate.py", f"Rastrigin{i}", f"pso_rastr{i}.py"])
