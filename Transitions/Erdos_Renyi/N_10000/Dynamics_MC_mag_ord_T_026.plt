set term postscript enhanced color eps dl 2.5
filenameoutput="Dynamics_MC_mag_ord_T_026.eps"
set output filenameoutput

reset
set size 1,1
set origin 0,0
set ylabel "{/=24 m}" rotate by 90 offset 0,0
set xlabel "{/=24 t}" offset 0,0
set key spacing 2 maxrows 5 width 6 at 9000, 0.8
set tics  font ",24"
set yrange[-0.05:1.05]
set xrange[-100:10500]
set ytics 0, 0.2, 1.0
set xtics 0, 2000, 10000

p "./MC_Erdos_Renyi_N_10000_k_3_mag_alpha_0.560_g_3_temp_0.260_tl_10000_tsampl_100_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.56}" w lp \
, "./MC_Erdos_Renyi_N_10000_k_3_mag_alpha_0.580_g_3_temp_0.260_tl_10000_tsampl_100_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.58}" w lp \
, "./MC_Erdos_Renyi_N_10000_k_3_mag_alpha_0.600_g_3_temp_0.260_tl_10000_tsampl_100_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.60}" w lp \
, "./MC_Erdos_Renyi_N_10000_k_3_mag_alpha_0.620_g_3_temp_0.260_tl_10000_tsampl_100_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.62}" w lp