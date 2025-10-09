set term postscript enhanced color eps dl 2.5
filenameoutput="Dynamics_MC_mag_rand_beta_250.eps"
set output filenameoutput

reset
set size 1,1
set multi
set origin 0,0
set ylabel "{/=24 m}" rotate by 90 offset 0,0
set xlabel "{/=24 t}" offset 0,0
set key spacing 2 maxrows 5 width 6 at 9000, 0.8
set tics  font ",24"

p "./Full_Random_Regular_Mag_deord_beta_2.5_nodes_10000_neigh_5_G_5_MCsteps_10000_sims_10_alpha_0.77_rand_parsed.txt" title "{/=24 {/Symbol a}=0.77}" w p \
, "./Full_Random_Regular_Mag_deord_beta_2.5_nodes_10000_neigh_5_G_5_MCsteps_10000_sims_10_alpha_0.8_rand_parsed.txt" title "{/=24 {/Symbol a}=0.80}" w p \
, "./Full_Random_Regular_Mag_deord_beta_2.5_nodes_10000_neigh_5_G_5_MCsteps_10000_sims_10_alpha_0.83_rand_parsed.txt" title "{/=24 {/Symbol a}=0.83}" w p \
, "./Full_Random_Regular_Mag_deord_beta_2.5_nodes_10000_neigh_5_G_5_MCsteps_10000_sims_10_alpha_0.87_rand_parsed.txt" title "{/=24 {/Symbol a}=0.87}" w p