# NSF_REU_Summer2019

This contains the NSGA-II swarm-algorithm evolutionary script that is used to evaluate and compare the multi-objective performance of such hybrid algorithms.
This work was supported by National Science Foundation Award (#1749635).

### Requirements
Requires Python 3.7 (or later) alongside numpy and scipy.
matplotlib may also be required if the visualization script is used.

### Running scripts
To run the model selection script, change the parameters within the script and use 
`python3 evaluate.py [Model] [output_filename]`
There are two models available, `Weibull` and `Covariate`. Including `-p` in the launch options can be used to evaluate pre-defined chromosome patterns in order to cross-validate data, and `-h` can be used to view objectives per-iteration in order to output to a file to plot histograms. NSGA-II parameters can be changed within the script. 

The `models.py` script contains objective functions for optimization, so a user may opt to modify or add their own swarm-applicable model. The swarm-hybrid algorithms will attempt to minimize whichever model it is given. This script is used by the `evaluate` script, and should not be ran on its own.

`python3 visualize.py [input_filename] [Model]` 
may also be used to plot various data (histograms, objectives, i.e. figures 5-9 in manuscript). Some slight editing of the script may be required to modify plot parameters. Included in the repository is the cross-validation data for the Covariate model testing, for both the histogram and test plotting.

## Algorithm Source
Algorithms involved were adapted from R implementations of Brownlee's *Clever Algorithms* ([Link](https://github.com/clever-algorithms/CleverAlgorithms)), most notably the NSGA-II-based overhead management, or *Nature-inspired metaheuristic algorithms* (Yang, Xin-She) in the case of the swarm algorithms involved.


## Performance
Accurate comparison between results depends on the performance of the computer running the script. The NSGA-II optimization is ran single-threaded in order to verify the integrity of all separate algorithm runs. As such, the script will take some time to run provided somewhat-large NSGA parameters. Algorithm tuning was performed using an Intel i9-9900K using stock clocks.

