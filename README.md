# A multi-threaded SPSA simulator in C

## Introduction

This program does multiple things

1. It simulates the chess version of the SPSA algorithm
invented by Joona Kiiski. This is a version of the SPSA algorithm
suitable for tuning parameters of chess engines. The simulation
depends on the Elo model together with a quadratic function `parameters-->elo`
which we will henceforth call the "(true) loss function".

2. It computes good (constant) SPSA parameters for a user requested final
"precision" (i.e. Elo below the optimum). This depends on the user
guessing the shape of a reasonable quadratic loss function. In other words the
user serves as an _oracle_. The computed parameters are
conservative however. They will also work for other reasonable loss
functions.

3. It evaluate theoretically the performance of the algorithm
with the given parameters and the true loss function. 
The theoretical characteristics can thus be compared to the simulated
ones (hint: there is perfect agreement).

Example output

```
$ ./spsa_sim --num_params 4

general options
===============
num_threads       =4
truncate          =-1
seed              =   1590606156
start_elo         =2.000000
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
elo average               =-0.215234 Elo
       fixed part         =-0.004958 Elo
       noise part         =-0.210276 Elo
50% percentile            =-0.180640 Elo
95% percentile            =-0.510440 Elo
success rate              =0.945835

sims
====

<output snipped>

sims=12852 success=0.9456[0.9396,0.9516] elo_avg=-0.214578 p50=0.498444 p95=0.950358

```


## The SPSA algorithm

## Choosing parameters and theoretical evaluation

### General principles

r <-> precision

rc^2 <-> speed at which the optimum is approached

