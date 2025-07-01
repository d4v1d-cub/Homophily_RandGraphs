import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

import scienceplots
plt.style.use(['science'])
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter  


def import_files(name,directory):
    dataframes = []
    files = [filename for filename in os.listdir(directory) if filename.startswith(name) and filename.endswith(".dat")]
    files.sort(key=lambda x: float(x.split('_')[-1].split('.')[0]+'.'+x.split('_')[-1].split('.')[1]))  # Sort based on the numerical part

    for filename in files:
            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as file:
                column_name = filename.split('_')[-1].split('.')[0]+'.'+filename.split('_')[-1].split('.')[1] # Use the file name as the column name
                data = [float(line.strip()) for line in file.readlines()]  # Read each line as a float number
                df = pd.DataFrame(data, columns=[column_name])
                dataframes.append(df)
    return dataframes


def mean_df(df,N):
    sims=len(df[0])/(N+1)
    mdfo=[]
    mdfdo=[]
    for i in range(len(df)):
        so=0
        sdo=0
        for j in range(int(sims)):
            if j<int(sims/2):
                so+=df[i].iloc[(N+1)*j:(N+1)*(j+1)].reset_index(drop=True)
            else:
                sdo+=df[i].iloc[(N+1)*j:(N+1)*(j+1)].reset_index(drop=True)
        mdfo.append(so/(sims/2))
        mdfdo.append(sdo/(sims/2))
    return mdfo,mdfdo


main_legend_fs=24
main_ylabel_fs=30
main_xlabel_fs=30
main_ticks_fs=28

minset_legend_fs=34
minset_ylabel_fs=36
minset_xlabel_fs=36
minset_ticks_fs=34

texts=30

#DATA - FIG. 1

filepath='data/sup_materials/fig_1/'
df_t_2 = pd.read_csv(filepath+'BP_RGfc_eq_c_3_T_0.2_delta_-10_maxiter_100000_ord_init_conf_1.0_G_2.txt', 
                 sep='\t', 
                 comment='#', 
                 names=['T', 'alpha', 'e', 'm(norm)', 'm(1)', 'free_ener', 'energy_dens', 'ρ+', 'convergence'])

df = pd.read_csv(filepath+'BP_RGfc_eq_c_3_T_0.5_delta_-10_maxiter_100000_ord_init_conf_1.0_G_2.txt', 
                 sep='\t', 
                 comment='#', 
                 names=['T', 'alpha', 'e', 'm(norm)', 'm(1)', 'free_ener','energy_dens', 'ρ+','convergence'])


#FIG. 1

print('Saving fig. 1 ...\n')

plt.figure(figsize=(12,8))
plt.plot(df['alpha'], df['ρ+'], 'r', linewidth=4,label='T=0.5')
plt.plot(df_t_2['alpha'], df_t_2['ρ+'], 'b', linewidth=4,label='T=0.2')
plt.xlabel(r'$\alpha$',fontsize=main_xlabel_fs)
plt.ylabel(r'$\rho +$',fontsize=main_ylabel_fs)
plt.gca().xaxis.set_major_locator(MaxNLocator(4))  
plt.gca().yaxis.set_major_locator(MaxNLocator(4))
plt.legend(fontsize=main_legend_fs, loc='lower right',handlelength=1)
plt.xticks(fontsize=main_ticks_fs)
plt.yticks(fontsize=main_ticks_fs)
plt.savefig('sm_figs/rho_+_g_2_T_0.2_0.5_c_3.pdf',dpi=800)


plt.figure(figsize=(12,8))
plt.plot(df['alpha'], df['energy_dens'], 'r', linewidth=4,label='T=0.5')
plt.plot(df_t_2['alpha'], df_t_2['energy_dens'], 'b', linewidth=4,label='T=0.2')


plt.xlabel(r'$\alpha$',fontsize=main_xlabel_fs)
plt.ylabel(r'$|e|$',fontsize=main_ylabel_fs)

plt.gca().xaxis.set_major_locator(MaxNLocator(4))  
plt.gca().yaxis.set_major_locator(MaxNLocator(4))
plt.legend(fontsize=main_legend_fs, loc='lower right',handlelength=1)
plt.xticks(fontsize=main_ticks_fs)
plt.yticks(fontsize=main_ticks_fs)

plt.savefig('sm_figs/e_dens_g_2_T_0.2_0.5_c_3.pdf',dpi=800)


#DATA - FIG.2


directory = "data/sup_materials/fig_2/"  

fmd=import_files("Full_Mag_deord_2",directory)
fmo=import_files("Full_Mag_ord_2",directory)

N=1000
mag_d_2=mean_df(fmd,N)
mag_o_2=mean_df(fmo,N)

with open(directory+"BP_RGfc_eq_N_1000_c_3_T_0.5_delta_-15_sampling_10_maxiter_1000_seed_1_G_2.txt", 'r') as file:
                frl= file.readlines()
                alphasx2=[float(frl[i].split()[1]) for i in range(len(frl))]
                magy2=[float(frl[i].split()[8]) for i in range(len(frl))]
                
                


with open(directory+"BP_RGfc_eq_N_1000_c_3_T_0.2_delta_-15_sampling_10_maxiter_1000_seed_1_G_2.txt", 'r') as file:
                frl= file.readlines()
                alphas=[float(frl[i].split()[1]) for i in range(len(frl))]
                mag=[float(frl[i].split()[8]) for i in range(len(frl))]
                
                
N=5000
fmd=import_files("Full_Mag_deord_5",directory)
fmo=import_files("Full_Mag_ord_5",directory)

mag_d=mean_df(fmd,N)
mag_o=mean_df(fmo,N)

y=[np.array(df[df.columns[0]])[-1] for df in mag_d[0]]
y2=[np.array(df[df.columns[0]])[-1] for df in mag_d[1]]

y3=[np.array(df[df.columns[0]])[-1] for df in mag_d_2[0]]
y4=[np.array(df[df.columns[0]])[-1] for df in mag_d_2[1]]


column_names = [float(df.columns[0]) for df in mag_d[0]]
column_names.sort()


#FIG. 2

print('Saving fig. 2 ...\n')


plt.figure(figsize=(12,8),dpi=1200)

plt.scatter(column_names,y,label='MC',marker='^',s=100,color='black', facecolor='blue',zorder=2)
plt.plot(alphas,mag,label='BP-T=0.2',color='blue',linewidth=3,zorder=1)

plt.scatter(column_names,y3,label='MC',marker='s',s=100,color='black', facecolor='red',zorder=2)
plt.plot(alphasx2,magy2,label='BP-T=0.5'  ,color='red',linewidth=3,zorder=1)


plt.xlabel(r'$\alpha$',fontsize=main_xlabel_fs)
plt.ylabel('m',fontsize=main_ylabel_fs)
plt.xticks(fontsize=main_ticks_fs)
plt.yticks(fontsize=main_ticks_fs)
plt.legend(fontsize=main_legend_fs)
plt.savefig('sm_figs/MagDvsalphas_1000and5000_T02T05_BP_G2.pdf',dpi=800)


#DATA - FIG. 3


print('Saving fig. 3 ...\n')

filepath='data/sup_materials/fig_3/'
data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_3.03_nodes_10000_neigh_4_G_4_MCsteps_10000_sims_10_alpha_0.65_ord_parsed.txt')

# Extract columns
x_beta_303_65 = data[:, 0]  # First column as x-axis
y_beta_303_65 = data[:, 1]  # Second column as y-axis


data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_3.03_nodes_10000_neigh_4_G_4_MCsteps_10000_sims_10_alpha_0.7_ord_parsed.txt')

x_beta_303_7 = data[:, 0]  
y_beta_303_7 = data[:, 1]  

plt.figure(figsize=(12, 8))

scatter1 = plt.scatter(x_beta_303_65, y_beta_303_65, 
                      marker='s', facecolors='none', 
                      edgecolors='lightcoral', s=100, linewidth=2,
                      label=r'$\alpha=0.65$')

scatter2 = plt.scatter(x_beta_303_7, y_beta_303_7, 
                     marker='^', facecolors='none', 
                     edgecolors='royalblue', s=100, linewidth=2,
                     label=r'$\alpha=0.7$')

plt.xlabel('t', fontsize=minset_xlabel_fs)
plt.ylabel('m', fontsize=minset_ylabel_fs)

plt.gca().xaxis.set_major_locator(MaxNLocator(5))  
plt.gca().yaxis.set_major_locator(MaxNLocator(5))
plt.ticklabel_format(axis='x', style='plain')

ax = plt.gca()

plt.xticks(fontsize=minset_ticks_fs)
plt.yticks(fontsize=minset_ticks_fs)

plt.legend(fontsize=minset_legend_fs)

plt.tight_layout()

plt.savefig('sm_figs/mag_vs_t_MC_303_ord.pdf',dpi=800)

data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_5.26_nodes_10000_neigh_4_G_4_MCsteps_100000_sims_10_alpha_0.54_ord_parsed.txt')

x_beta_526_54 = data[:, 0]
y_beta_526_54 = data[:, 1]


data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_5.26_nodes_10000_neigh_4_G_4_MCsteps_100000_sims_10_alpha_0.59_ord_parsed.txt')

x_beta_526_59 = data[:, 0] 
y_beta_526_59 = data[:, 1] 



plt.figure(figsize=(12, 8))

scatter1 = plt.scatter(x_beta_526_54, y_beta_526_54, 
                      marker='s', facecolors='none', 
                      edgecolors='lightcoral', s=100, linewidth=2,
                      label=r'$\alpha=0.54$')

scatter2 = plt.scatter(x_beta_526_59, y_beta_526_59, 
                     marker='^', facecolors='none', 
                     edgecolors='royalblue', s=100, linewidth=2,
                     label=r'$\alpha=0.59$')

plt.gca().xaxis.set_major_locator(MaxNLocator(5))
plt.gca().yaxis.set_major_locator(MaxNLocator(5))

plt.ticklabel_format(axis='x', style='plain')
ax = plt.gca()

plt.xlabel('t', fontsize=minset_xlabel_fs)
plt.ylabel('m', fontsize=minset_ylabel_fs)

plt.xticks(fontsize=minset_ticks_fs)
plt.yticks(fontsize=minset_ticks_fs)

plt.legend(fontsize=minset_legend_fs)

plt.tight_layout()

plt.savefig('sm_figs/mag_vs_t_MC_526_ord.pdf',dpi=800)


data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_3.03_nodes_10000_neigh_4_G_4_MCsteps_10000_sims_10_alpha_0.67_rand_parsed.txt')

x_beta_303_67 = data[:, 0]  
y_beta_303_67 = data[:, 1]  


data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_3.03_nodes_10000_neigh_4_G_4_MCsteps_10000_sims_10_alpha_0.7_rand_parsed.txt')

x_beta_303_7 = data[:, 0]  
y_beta_303_7 = data[:, 1]  

data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_3.03_nodes_10000_neigh_4_G_4_MCsteps_10000_sims_10_alpha_0.73_rand_parsed.txt')

x_beta_303_73 = data[:, 0]  
y_beta_303_73 = data[:, 1]  


data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_3.03_nodes_10000_neigh_4_G_4_MCsteps_10000_sims_10_alpha_0.78_rand_parsed.txt')

x_beta_303_78 = data[:, 0]
y_beta_303_78 = data[:, 1]


plt.figure(figsize=(12, 8))

scatter1 = plt.scatter(x_beta_303_67, y_beta_303_67, 
                      marker='s', facecolors='none', 
                      edgecolors='lightcoral', s=100, linewidth=2,
                      label=r'$\alpha=0.67$')

scatter3 = plt.scatter(x_beta_303_7, y_beta_303_7, 
                      marker='s', facecolors='none', 
                      edgecolors='peachpuff', s=100, linewidth=2,
                      label=r'$\alpha=0.70$')

scatter4 = plt.scatter(x_beta_303_73, y_beta_303_73, 
                      marker='s', facecolors='none', 
                      edgecolors='lightgreen', s=100, linewidth=2,
                      label=r'$\alpha=0.73$')


scatter2 = plt.scatter(x_beta_303_78, y_beta_303_78, 
                     marker='^', facecolors='none', 
                     edgecolors='royalblue', s=100, linewidth=2,
                     label=r'$\alpha=0.78$')

plt.gca().xaxis.set_major_locator(MaxNLocator(5))  # 5 ticks on x-axis
plt.gca().yaxis.set_major_locator(MaxNLocator(5))

plt.ticklabel_format(axis='x', style='plain')
ax = plt.gca()

plt.xlabel('t', fontsize=minset_xlabel_fs)
plt.ylabel('m', fontsize=minset_ylabel_fs)

plt.xticks(fontsize=minset_ticks_fs)
plt.yticks(fontsize=minset_ticks_fs)

plt.legend(fontsize=minset_legend_fs)

plt.tight_layout()

plt.savefig('sm_figs/mag_vs_t_MC_303_rand.pdf',dpi=800)



data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_5.26_nodes_10000_neigh_4_G_4_MCsteps_100000_sims_10_alpha_0.54_rand_parsed.txt')

x_beta_526_54 = data[:, 0]  
y_beta_526_54 = data[:, 1]  

data = np.loadtxt(filepath+'Full_Random_Regular_Mag_deord_beta_5.26_nodes_10000_neigh_4_G_4_MCsteps_100000_sims_10_alpha_0.59_rand_parsed.txt')

x_beta_526_59 = data[:, 0]  
y_beta_526_59 = data[:, 1]  

plt.figure(figsize=(12, 8))

scatter1 = plt.scatter(x_beta_526_54, y_beta_526_54, 
                      marker='s', facecolors='none', 
                      edgecolors='lightcoral', s=100, linewidth=2,
                      label=r'$\alpha=0.54$')

scatter2 = plt.scatter(x_beta_526_59, y_beta_526_59, 
                     marker='^', facecolors='none', 
                     edgecolors='royalblue', s=100, linewidth=2,
                     label=r'$\alpha=0.59$')

plt.gca().xaxis.set_major_locator(MaxNLocator(5))  
plt.gca().yaxis.set_major_locator(MaxNLocator(5))

plt.ticklabel_format(axis='x', style='plain')
ax = plt.gca()

plt.xlabel('t', fontsize=minset_xlabel_fs)
plt.ylabel('m', fontsize=minset_ylabel_fs)

plt.xticks(fontsize=minset_ticks_fs)
plt.yticks(fontsize=minset_ticks_fs)

plt.legend(fontsize=minset_legend_fs)

plt.tight_layout()

plt.savefig('sm_figs/mag_vs_t_MC_526_rand.pdf',dpi=800)


#DATA - FIG. 5

print('Saving fig. 5 ...\n')

file_path='data/paper/paper_fig_3/inset/'
file = file_path+'BP_ER_popdyn_eq_c_3_T_0.200_delta_1e-6_maxiter_1000_init_p_0.51_G_3_deltamag_0.01_popsize_100000_smoothsize_100_measuresize_1000_seed_1.txt'



data = np.loadtxt(file)

f_x = data[:, 1]  
f_y = data[:, 4]  

file = file_path+'BP_ER_popdyn_eq_c_3_T_0.200_delta_1e-6_maxiter_1000_init_p_1.0_G_3_deltamag_0.01_popsize_100000_smoothsize_100_measuresize_1000_seed_1.txt'



data = np.loadtxt(file)

f_x_ord= data[:, 1] 
f_y_ord = data[:, 4]


plt.figure(figsize=(12, 8), dpi=600, facecolor='none')

plt.scatter(f_x,f_y, 
            color='red', label='BP-rand', marker='s', edgecolor='red', facecolor='none',linewidths=1,s=50)  

plt.scatter(f_x_ord, f_y_ord, 
            color='blue', label='BP-ord', marker='^', edgecolor='blue', facecolor='none',linewidths=1,s=50)  



ax = plt.gca()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.legend(fontsize=minset_legend_fs, markerscale=2)
plt.ylabel('f', fontsize=minset_ylabel_fs)
plt.xlabel(r'$\alpha$', fontsize=minset_xlabel_fs)
plt.tick_params(axis='both', labelsize=minset_ticks_fs)
plt.text(0.9, -6, 'T = 0.2', rotation=0, fontsize=texts, color='black')
plt.tight_layout()

plt.savefig('sm_figs/er_free_energy_T_0.2' + '.pdf', dpi=800)

file_path='data/sup_materials/fig_5/'
file = file_path+'BP_ER_popdyn_eq_c_3_T_0.500_delta_1e-6_maxiter_1000_init_p_0.51_G_3_deltamag_0.01_popsize_100000_smoothsize_100_measuresize_1000_seed_1.txt'



data = np.loadtxt(file)

f_x = data[:, 1]  
f_y = data[:, 4]  
file = file_path+'BP_ER_popdyn_eq_c_3_T_0.500_delta_1e-6_maxiter_1000_init_p_1.0_G_3_deltamag_0.01_popsize_100000_smoothsize_100_measuresize_1000_seed_1.txt'



data = np.loadtxt(file)

f_x_ord= data[:, 1]  
f_y_ord = data[:, 4] 

plt.figure(figsize=(12, 8), dpi=600, facecolor='none')

plt.scatter(f_x,f_y, 
            color='red', label='BP-rand', marker='s', edgecolor='red', facecolor='none',linewidths=1,s=50)  

plt.scatter(f_x_ord, f_y_ord, 
            color='blue', label='BP-ord', marker='^', edgecolor='blue', facecolor='none',linewidths=1,s=50)  


ax = plt.gca()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.legend(fontsize=minset_legend_fs, markerscale=2)
plt.ylabel('f', fontsize=minset_xlabel_fs)
plt.xlabel(r'$\alpha$', fontsize=minset_xlabel_fs)
plt.tick_params(axis='both', labelsize=minset_ticks_fs)
plt.text(0.6, -3.2, 'T = 0.5', rotation=0, fontsize=texts, color='black')

plt.savefig('sm_figs/er_free_energy_T_0.5' + '.pdf', dpi=800)



file_path='data/paper/paper_fig_1/inset/'
file = file_path+'BP_RGfc_allalpha_eq_c_5_delta_-10_maxiter_100000_hom_init_p_0.51_G_5_alphas_0.5_0.005_1.0_temps_0.05_0.01_1.2.txt'



data = np.genfromtxt(file, delimiter='\t', skip_header=1, dtype=None, encoding=None)

alpha_rand = data['f1']  
f_rand = data['f5']


file = file_path+'BP_RGfc_allalpha_eq_c_5_delta_-10_maxiter_100000_hom_init_p_1.0_G_5_alphas_0.5_0.005_1.0_temps_0.05_0.01_1.2.txt'



data = np.genfromtxt(file, delimiter='\t', skip_header=1, dtype=None, encoding=None)

alpha_ord = data['f1']  
f_ord = data['f5']    



plt.figure(figsize=(12, 8), dpi=600, facecolor='none')


plt.scatter(alpha_rand[3031:3131], f_rand[3031:3131], 
            color='red', label='BP-rand', marker='s', edgecolor='red', facecolor='none',linewidths=1,s=50) 
plt.scatter(alpha_ord[3031:3131], f_ord[3031:3131], 
            color='blue', label='BP-ord', marker='^', edgecolor='blue', facecolor='none',linewidths=1,s=50)  
ax = plt.gca()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

plt.legend(fontsize=minset_legend_fs, markerscale=2)
plt.ylabel('f', fontsize=minset_ylabel_fs)
plt.xlabel(r'$\alpha$', fontsize=minset_xlabel_fs)
plt.tick_params(axis='both', labelsize=minset_ticks_fs)
plt.text(0.6, -2.3, 'T = 0.35', rotation=0, fontsize=texts, color='black')

plt.savefig('sm_figs/free_energy_T_0.35' + '.pdf', dpi=800)

file_path='data/sup_materials/fig_5/'

file = file_path+'BP_RGfc_allalpha_eq_c_5_delta_-10_maxiter_100000_hom_init_p_0.51_G_5_alphas_0.5_0.005_1.0_temps_0.15_0.01_0.15.txt'



data = np.genfromtxt(file, delimiter='\t', skip_header=1, dtype=None, encoding=None)

alpha_rand = data['f1'] 
f_rand = data['f5']

file = file_path+'BP_RGfc_allalpha_eq_c_5_delta_-10_maxiter_100000_hom_init_p_1.0_G_5_alphas_0.5_0.005_1.0_temps_0.15_0.01_0.15.txt'



data = np.genfromtxt(file, delimiter='\t', skip_header=1, dtype=None, encoding=None)

alpha_ord = data['f1'] 
f_ord = data['f5']    

plt.figure(figsize=(12, 8), dpi=600, facecolor='none')

plt.scatter(alpha_rand, f_rand, 
            color='red', label='BP-rand', marker='s', edgecolor='red', facecolor='none',linewidths=1,s=50)  
plt.scatter(alpha_ord, f_ord, 
            color='blue', label='BP-ord', marker='^', edgecolor='blue', facecolor='none',linewidths=1,s=50)  
ax = plt.gca()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

plt.legend(fontsize=minset_legend_fs, markerscale=2)
plt.ylabel('f', fontsize=minset_ylabel_fs)
plt.xlabel(r'$\alpha$', fontsize=minset_xlabel_fs)
plt.tick_params(axis='both', labelsize=minset_ticks_fs)
plt.text(0.6, -2, 'T = 0.15', rotation=0, fontsize=texts, color='black')


plt.savefig('sm_figs/free_energy_T_0.15' + '.pdf', dpi=800)
