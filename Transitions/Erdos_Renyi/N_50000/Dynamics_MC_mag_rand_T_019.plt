set term postscript enhanced color eps dl 2.5
filenameoutput="Dynamics_MC_mag_rand_T_019.eps"
set output filenameoutput

reset
set size 1,1
set origin 0,0
set ylabel "{/=24 m}" rotate by 90 offset 0,0
set xlabel "{/=24 t}" offset 0,0
set key spacing 2 maxrows 5 width 6 at 90000, 0.8
set tics  font ",24"
set yrange[-0.05:1.05]
set xrange[-100:105000]
set ytics 0, 0.2, 1.0
set xtics 0, 20000, 100000

p "./MC_Erdos_Renyi_N_50000_k_3_mag_alpha_0.520_g_3_temp_0.190_tl_100000_tsampl_1000_sgraph0_1_shist0_1_nhist_1_init_rand.txt" title "{/=24 {/Symbol a}=0.52}" w lp \
, "./MC_Erdos_Renyi_N_50000_k_3_mag_alpha_0.540_g_3_temp_0.190_tl_100000_tsampl_1000_sgraph0_1_shist0_1_nhist_1_init_rand.txt" title "{/=24 {/Symbol a}=0.54}" w lp \
, "./MC_Erdos_Renyi_N_50000_k_3_mag_alpha_0.560_g_3_temp_0.190_tl_100000_tsampl_1000_sgraph0_1_shist0_1_nhist_1_init_rand.txt" title "{/=24 {/Symbol a}=0.56}" w lp \
, "./MC_Erdos_Renyi_N_50000_k_3_mag_alpha_0.580_g_3_temp_0.190_tl_100000_tsampl_1000_sgraph0_1_shist0_1_nhist_1_init_rand.txt" title "{/=24 {/Symbol a}=0.58}" w lp