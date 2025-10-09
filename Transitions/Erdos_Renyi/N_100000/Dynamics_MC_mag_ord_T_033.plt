set term postscript enhanced color eps dl 2.5
filenameoutput="Dynamics_MC_mag_ord_T_033.eps"
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

p "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.590_g_3_temp_0.330_tl_1000_tsampl_10_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.59}" w lp \
, "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.610_g_3_temp_0.330_tl_1000_tsampl_10_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.61}" w lp \
, "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.620_g_3_temp_0.330_tl_10000_tsampl_100_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.62}" w lp \
, "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.630_g_3_temp_0.330_tl_10000_tsampl_100_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.63}" w lp \
, "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.640_g_3_temp_0.330_tl_10000_tsampl_100_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.64}" w lp \
, "./MC_Erdos_Renyi_N_100000_k_3_mag_alpha_0.650_g_3_temp_0.330_tl_1000_tsampl_10_sgraph0_1_shist0_1_nhist_1_init_ord.txt" title "{/=24 {/Symbol a}=0.65}" w lp