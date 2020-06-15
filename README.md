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
seed              =   1592246463274356340
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
sims=91 success=0.9451[0.8734,1.0167] elo_avg=-0.234242 p50=0.505495 p95=0.945055
sims=187 success=0.9519[0.9049,0.9988] elo_avg=-0.214012 p50=0.556150 p95=0.957219
sims=286 success=0.9615[0.9274,0.9957] elo_avg=-0.203717 p50=0.559441 p95=0.965035

...
sims=8786 success=0.9456[0.9383,0.9529] elo_avg=-0.214955 p50=0.500228 p95=0.949237
sims=8885 success=0.9456[0.9384,0.9529] elo_avg=-0.214848 p50=0.500619 p95=0.949240
sims=8971 success=0.9458[0.9387,0.9530] elo_avg=-0.214844 p50=0.499944 p95=0.949392
sims=9068 success=0.9459[0.9387,0.9530] elo_avg=-0.214762 p50=0.500221 p95=0.949382
sims=9163 success=0.9461[0.9390,0.9532] elo_avg=-0.214542 p50=0.500491 p95=0.949689
...

```



