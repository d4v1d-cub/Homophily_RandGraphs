set term postscript enhanced color eps dl 2.5
filenameoutput="Dynamics_MC_mag_rand_beta_213.eps"
set output filenameoutput

reset
set size 1,1
set multi
set origin 0,0
set ylabel "{/=24 m}" rotate by 90 offset 0,0
set xlabel "{/=24 t}" offset 0,0
set key spacing 2 maxrows 5 width 6 at 90000, 0.8
set tics  font ",24"

p "./Full_Random_Regular_Mag_deord_beta_2.13_nodes_1000_neigh_5_G_5_MCsteps_100000_sims_10_alpha_0.86_rand_parsed.txt" title "{/=24 {/Symbol a}=0.86}" w p \
, "./Full_Random_Regular_Mag_deord_beta_2.13_nodes_1000_neigh_5_G_5_MCsteps_100000_sims_10_alpha_0.89_rand_parsed.txt" title "{/=24 {/Symbol a}=0.89}" w p \
, "./Full_Random_Regular_Mag_deord_beta_2.13_nodes_1000_neigh_5_G_5_MCsteps_100000_sims_10_alpha_0.92_rand_parsed.txt" title "{/=24 {/Symbol a}=0.92}" w p \
, "./Full_Random_Regular_Mag_deord_beta_2.13_nodes_1000_neigh_5_G_5_MCsteps_100000_sims_10_alpha_0.96_rand_parsed.txt" title "{/=24 {/Symbol a}=0.96}" w p