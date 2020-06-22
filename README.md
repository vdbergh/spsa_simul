# A multi-threaded SPSA simulator in C

## Introduction

This program does multiple things

1. It simulates the chess version of the SPSA algorithm
invented by Joona Kiiski. This is a version of the SPSA algorithm
suitable for tuning parameters of chess engines. The simulation
depends on the Elo model together with a quadratic function `parameters-->elo`
which we will henceforth call the "true loss function".

2. It computes good (constant) SPSA parameters for a user requested
final "precision" (i.e. Elo below the optimum). This depends on the
user guessing the shape of a reasonable quadratic loss function,
called henceforth the "guessed loss function". In other words the user
serves as an _oracle_. The computed parameters are conservative
however. They will also work for other reasonable loss functions.

3. It evaluate theoretically the performance of the algorithm with the
given parameters and the true loss function.  The theoretical
characteristics can thus be compared to the simulated ones (hint:
there is perfect agreement).

For the theoretical background see

https://github.com/vdbergh/spsa_simul/tree/master/doc/theoretical_basis.pdf

## Example output

```
$ ./spsa_sim --num_params 4

general options
===============
num_threads       =4
truncate          =-1
true_start_elo    =2.000000
heuristic         =1
seed              =   1592808355013632880
quiet             =0

guessed loss function
=====================
+---------------+-----------------+-----------------+---------------+
|      elos     |      minima     |      optima     |     maxima    |
+---------------+-----------------+-----------------+---------------+
|        2.00   |         0.00    |       100.00    |       200.00  |
|        2.00   |         0.00    |       100.00    |       200.00  |
|        2.00   |         0.00    |       100.00    |       200.00  |
|        2.00   |         0.00    |       100.00    |       200.00  |
+---------------+-----------------+-----------------+---------------+

true loss function
==================
+---------------+-----------------+-----------------+---------------+
|      elos     |      minima     |      optima     |     maxima    |
+---------------+-----------------+-----------------+---------------+
|        2.00   |         0.00    |       100.00    |       200.00  |
|        2.00   |         0.00    |       100.00    |       200.00  |
|        2.00   |         0.00    |       100.00    |       200.00  |
|        2.00   |         0.00    |       100.00    |       200.00  |
+---------------+-----------------+-----------------+---------------+

spsa data
=========
~~~~~~~design~~~~~~~ 
num_params        =4
confidence        =0.950000
draw_ratio        =0.610000
precision         =0.500000
c_ratio           =0.166667
lambda_ratio      =3.000000
start_elo         =2.000000
bounds            =0
~~~~~~~computed~~~~~~~ 
r                 =0.003111
c                 =33.333333  33.333333  33.333333  33.333333  
num_games         =376868

starting point sims
===================
starting point      =150.000000  150.000000  150.000000  150.000000  
true loss function  =-2.000000

theoretical characteristics (using the true loss function)
==========================================================
Elo half life (games)     = 43537.498565  43537.498565  43537.498565  43537.498565  
Elo average               = -0.215234 Elo (asymptotic: -0.210799, scale factor: 1.021039)
50% percentile            = -0.180640 Elo (asymptotic: -0.176897, scale factor: 1.021159)
95% percentile            = -0.510440 Elo (asymptotic: -0.500000, scale factor: 1.020879)
success rate              = 0.945835

sims
====
sims=91 success=0.9341[0.8560,1.0121] elo_avg=-0.193604 p50=0.549451 p95=0.934066

...

sims=3062 success=0.9376[0.9245,0.9507] elo_avg=-0.220409 p50=0.486937 p95=0.942521

```

## SPSA

For the version of SPSA we use see

https://github.com/zamar/spsa

Note that in contrast to loc. cit. we keep the hyperparameters `ck`
and `Rk:=ak/ck^2` constant and we denote them by `c` and `r`
respectively. When optimizing more than one parameter at a time, the
value of `c` is allowed to depend on the parameter but `r` is not.

The rationale for decreasing `ck` and `Rk` over time in the standard
SPSA algorithm is to achieve "convergence". However it can be easily
shown that in our setting, due to lack of resources, no reasonable
form of convergence is ever achievable. The best one can do is to plan
for ending up within a specified distance from the optimum (Elo wise)
with a specified probability. This is what what the current program
achieves, under some assumptions on the true loss function. The
specified distance is called `precision` (default `0.5 Elo`) and the
probability is called `confidence` (default `95%`).

## Parameters

```
$ ./spsa_sim -h
spsa_sim [-h] [--num_params NUM_PARAMS] [--confidence CONFIDENCE] [--draw_ratio DRAW_RATIO] [--seed SEED] [--truncate TRUNCATE] [--bounds] [--precision PRECISION] [--c_ratio C_RATIO] [--lambda_ratio LAMBDA_RATIO] [--est_elos EST_ELOS1,...] [--true_elos TRUE_ELOS1,...] [--minima MINIMA1,..] [--optima OPTIMA1,..] [--maxima MAXIMA1,...] [--start_elo START_ELO] [--quiet] [--threads THREADS]
```

### Lists

Some arguments have to be specified for every parameter. This is
done by supplying a comma (or colon) separated list, with no spaces.

For example:

```
--true_elos 1,2,3
```

By convention the last entry in such a list is repeated if the list is
shorter than the number of parameters.

### Specification of the loss functions

Both the true loss function and the guessed loss function are assumed
to be defined on a hypercube with edges `[minimum_k,maximum_k]` where
`_k` refers to the `k`'th parameter. Moreover they are also assumed
to achieve their optimum in the same point. Note: this is a dubious
assumption which needs to be changed.

*Arguments:*

```
--minima <list>
--maxima <list>
--optima <list>
```
The defaults are respectively `0,200,100`.

The actual loss functions are fixed by specifying Elo lists

```
--est_elos <list>
--true_elos <list?
```
The defaults is  `2` for both loss functions.

If `elo_k` is an entry in such a list then this means an Elo loss of `elo_k`
if the parameter is at a distance `(maximum_k-minimum_k)/2` from `optimum_k`. 

### Other parameters


```
--num_params NUM_PARAMS
```
The number of parameters being tuned (default: `1`).

```
--precision PRECISION
```
Requested Elo distance from the optimum (default: `0.5`).

```
--confidence CONFIDENCE
```
The probability of achieving the requested precision (default `0.95`).

```
--draw_ratio DRAW_RATIO
```
The draw_ratio (default `0.61`).

```
--seed SEED
```
The prng seed. This is a 64 bit number, expressed in decimal notation. If this
argument is not supplied then the prng is seeded using the time.

```
--truncate TRUNCATE
```
Stop simulating after this many sims. By default the simulator runs forever.

```
--bounds
```
When this flag is present the loss function will not be evaluated
outside the specified hypercube.  In that case the results of the
simulations will deviate slightly from the theoretical predictions,
which assume that the loss function can be evaluated everywhere.


```
--c_ratio C_RATIO
```
Set `c` (for a given parameter) equal to `(maximum-minum)*C_RATIO` (default: 1/6).

```
--lambda_ratio LAMBDA_RATIO
```
Rougly speaking `lambda` is the number of games it takes to divide the
Elo distance to the optimum by `e^2=7.389`. The total number of games
will be lambda*LAMBDA_RATIO (default: 3).

```
--start_elo START_ELO
```
Select a starting point for the simulation where the guessed loss
function takes the value START_ELO (default: 2).

```
--threads THREADS
```
Simultaneous simulations. The default is given by the number of cores
of on the system, as returned by the `nproc` command.

```
--quiet
```
When this flag is present very little information is printed.

## Output

The most essential output is given by the values of `r`, `c` and the
number of games.

_Example_

```
$ ./spsa_sim --num_params 4 --quiet --truncate 0
num_params        =4
start_elo         =2.00
num_games         =376868
r                 =0.003111
c                 =33.333333  33.333333  33.333333  33.333333  
```

