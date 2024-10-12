#!/bin/bash

exponent=-10

alpha_min=1.0
d_alpha=0.05
alpha_max=1.0

temp_min=0.05
d_temp=0.005
temp_max=1.2

gmin=1
gmax=4

cond_init=ord
# param=0.51
param=1

for g in $(seq $gmin $gmax)
do

c=3

julia map_states.jl $c $g $exponent $alpha_min $d_alpha $alpha_max $temp_min $d_temp $temp_max $cond_init $param

done