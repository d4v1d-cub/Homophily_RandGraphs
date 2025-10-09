set term postscript enhanced color eps dl 2.5
filenameoutput="Dynamics_MC_mag_rand_T_040.eps"
set output filenameoutput

reset
set size 1,1
set origin 0,0
set ylabel "{/=24 m}" rotate by 90 offset 0,0
set xlabel "{/=24 t}" offset 0,0
set key spacing 2 maxrows 5 width 6 at 900, 0.8
set tics  font ",24"
set yrange[-0.05:1.05]
set xrange[-100:1050]
set ytics 0, 0.2, 1.0
set xtics 0, 200, 1000

p "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.730_g_3_temp_0.400_tl_1000_tsampl_10_sgraph0_1_shist0_1_nhist_1_init_rand.txt" title "{/=24 {/Symbol a}=0.73}" w lp \
, "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.750_g_3_temp_0.400_tl_1000_tsampl_10_sgraph0_1_shist0_1_nhist_1_init_rand.txt" title "{/=24 {/Symbol a}=0.75}" w lp \
, "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.770_g_3_temp_0.400_tl_1000_tsampl_10_sgraph0_1_shist0_1_nhist_1_init_rand.txt" title "{/=24 {/Symbol a}=0.77}" w lp \
, "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.790_g_3_temp_0.400_tl_1000_tsampl_10_sgraph0_1_shist0_1_nhist_1_init_rand.txt" title "{/=24 {/Symbol a}=0.79}" w lp \
, "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.810_g_3_temp_0.400_tl_1000_tsampl_10_sgraph0_1_shist0_1_nhist_1_init_rand.txt" title "{/=24 {/Symbol a}=0.81}" w lp