import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import seaborn as sns
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.image as mpimg
from scipy.optimize import curve_fit

from matplotlib.ticker import FormatStrFormatter  
from matplotlib.ticker import ScalarFormatter, MaxNLocator


import re
def mean_dic_err(dic):
    l=[]
    error=[]
    beta=[]
    for key in dic.keys():
        mu=sum(dic[key])/len(dic[key])
        err=(dic[key][0]-dic[key][1])/2
        l.append(mu)
        error.append(err)
        beta.append(1/float(key))
    return l, beta, error


def gap_comp_int(array1,array2):
    delta=np.array(array2)-np.array(array1)
    if sum(delta)==0:
        i1=0
        i2=len(array1)
    else:
        i1=np.nonzero(delta)[0][0]
        i2=np.where(array2==1)[0][0]
    return sum(delta[i1:i2])/len(delta[i1:i2]) 

def calculate_area(x_list, y_list):
    
    area = 0.0
    for i in range(len(x_list) - 1):
        dx = x_list[i+1] - x_list[i]
        avg_height = ( (1 - y_list[i]) + (1 - y_list[i+1]) ) / 2
        area += dx * avg_height
    return area

def gap_int(array1,array2):
    delta=np.array(array2)-np.array(array1)
    if sum(delta)==0:
        i1=0
        i2=len(array1)
    else:
        i1=np.nonzero(delta)[0][0]
        i2=np.where(array2==1)[0][0]
    return i1,i2        

import scienceplots
plt.style.use(['science'])


main_legend_fs=24
main_ylabel_fs=30
main_xlabel_fs=30
main_ticks_fs=28

inset_legend_fs=22
inset_ylabel_fs=27
inset_xlabel_fs=27
inset_ticks_fs=24


#LOADING DATA


file_path='data/paper/paper_fig_1/main/'

file = file_path+'Phase_Diagram_Homophily_RGfc_c_5_G_5.txt'



# Read the data from the file
data = np.loadtxt(file)

# Separate the data into two arrays
xg5 = data[:, 0]  # First column
yg5 = data[:, 1]  # Second column

x=np.arange(0.01,0.55,0.01)


mc_dic_ord={"2.00":[0.9,0.86],"2.13":[0.85,0.81],"2.50":[0.75,0.71],"3.03":[0.66,0.62],"3.85":[0.59,0.55],"5.26":[0.57,0.53]}
mc_dic_rand={"2.00":[1.00,0.96],"2.13":[0.96,0.92],"2.50":[0.87,0.83],"3.03":[0.79,0.75],"3.85":[0.69,0.65],"5.26":[0.57,0.53]}

mc_ord=mean_dic_err(mc_dic_ord)
mc_rand=mean_dic_err(mc_dic_rand)

file_pattern = file_path+'Phase_Diagram_Homophily_numerical_c_5*.txt'
files = sorted(glob.glob(file_pattern))
files


all_y=[]

for file in files:
    data= np.loadtxt(file)
    y=data[:,1]
    all_y.append(y)


from mpl_toolkits.axes_grid1 import make_axes_locatable



x=np.arange(0.01,0.55,0.01)

mc_dic_ord={"2.00":[0.9,0.86],"2.13":[0.85,0.81],"2.50":[0.75,0.71],"3.03":[0.66,0.62],"3.85":[0.59,0.55],"5.26":[0.57,0.53]}
mc_dic_rand={"2.00":[1.00,0.96],"2.13":[0.96,0.92],"2.50":[0.87,0.83],"3.03":[0.79,0.75],"3.85":[0.69,0.65],"5.26":[0.57,0.53]}

mc_ord=mean_dic_err(mc_dic_ord)
mc_rand=mean_dic_err(mc_dic_rand)


file_pattern = file_path + 'Phase_Diagram_Homophily_numerical_c_5_delta_-10_maxiter_100000_*.txt'
files = sorted(glob.glob(file_pattern), key=lambda x: float(x.split('_')[-3]))

all_yG5 = []

# Read each file and extract the second column
for file in files:
    data = np.loadtxt(file)
    y = data[:, 1]
    all_yG5.append(y)


file_path='data/paper/paper_fig_1/inset/'

file = file_path+'BP_RGfc_allalpha_eq_c_5_delta_-10_maxiter_100000_hom_init_p_0.51_G_5_alphas_0.5_0.005_1.0_temps_0.05_0.01_1.2.txt'



data = np.genfromtxt(file, delimiter='\t', skip_header=1, dtype=None, encoding=None)

# Extract the second and third columns
alpha_rand = data['f1']  
f_rand = data['f5']

file = file_path+'BP_RGfc_allalpha_eq_c_5_delta_-10_maxiter_100000_hom_init_p_1.0_G_5_alphas_0.5_0.005_1.0_temps_0.05_0.01_1.2.txt'



data = np.genfromtxt(file, delimiter='\t', skip_header=1, dtype=None, encoding=None)

# Extract the second and third columns
alpha_ord = data['f1']  
f_ord = data['f5']    


#FIG. 1

plt.figure(figsize=(12, 8), dpi=600)


main_ax = plt.gca()  

cmap = plt.get_cmap('coolwarm_r', len(all_yG5))
for i in range(len(all_yG5)):
    if i==0 or i==len(all_yG5)-1:
        main_ax.plot(x, all_yG5[i], color=cmap(i), linewidth=3)
    else:
        main_ax.plot(x, all_yG5[i], color=cmap(i), linewidth=3)

main_ax.plot(xg5, yg5, label='Prediction', color='black', linestyle='--', linewidth=3)
main_ax.errorbar(mc_rand[1], mc_rand[0], yerr=mc_rand[2], fmt='s', markersize=12,
             markerfacecolor='red', color='black', label='MC-rand',
             capsize=5, elinewidth=2)
main_ax.errorbar(mc_ord[1], mc_ord[0], yerr=mc_ord[2], fmt='^', markersize=10,
             markerfacecolor='blue', color='black', label='MC-ord',
             capsize=5, elinewidth=2)

# Main plot formatting
main_ax.set_ylabel(r'$\alpha_c$', fontsize=30)
main_ax.set_xlabel("T", fontsize=30, labelpad=-8)
main_ax.tick_params(axis='both', labelsize=28)
main_ax.set_ylim(0.5, 1.025)
main_ax.legend(fontsize=24, loc='lower right')


# Add text annotations
main_ax.text(0.23, 0.70, 'Consensus is easy', rotation=37, fontsize=38, color='black')
main_ax.text(0.27, 0.62, 'Consensus is hard', rotation=40, fontsize=35, color='black')
main_ax.text(0.28, 0.52, 'Consensus is impossible', rotation=40, fontsize=35, color='black')

inset_ax = inset_axes(main_ax, 
                     width="60%", 
                     height="60%",
                     loc='upper left',
                     bbox_to_anchor=(0.08, 0.3, 0.7, 0.7),
                     bbox_transform=main_ax.transAxes)

# Plot inset content
inset_ax.scatter(alpha_rand[3031:3131], f_rand[3031:3131], 
                color='red', label='BP-rand', marker='s', edgecolor='red', 
                facecolor='none', linewidths=1, s=20)
inset_ax.scatter(alpha_ord[3031:3131], f_ord[3031:3131], 
                color='blue', label='BP-ord', marker='^', edgecolor='blue', 
                facecolor='none', linewidths=1, s=20)

# Format inset
inset_ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
inset_ax.legend(fontsize=22, markerscale=1.5, handletextpad=0.3)
inset_ax.set_ylabel('f', fontsize=27, labelpad=-7)
inset_ax.set_xlabel(r'$\alpha$', fontsize=27, labelpad=-3)
inset_ax.tick_params(axis='both', labelsize=24)
inset_ax.text(0.85, -1.9, 'T = 0.35', rotation=0, fontsize=24, color='black')


divider = make_axes_locatable(main_ax)
cax = divider.append_axes("bottom", size="5%", pad=0.8)
norm = mpl.colors.Normalize(vmin=0, vmax=len(all_yG5) - 1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
tick_indices = np.linspace(0, len(all_yG5) - 1, 5).astype(int)
cbar.ax.set_ylabel('p', fontsize=24, rotation=0, labelpad=20)
cbar.set_ticks(tick_indices)
cbar.set_ticklabels([files[i].split('_')[-3] for i in tick_indices])
cbar.ax.tick_params(labelsize=24)


plt.savefig('paper_figs/fig_1_phase_diagram_free_e_T_0.35' + '.jpg', dpi=800)




#LOADING DATA FIG 2
all_xs_c4 = []
all_ys_c4 = []

x_data=[1,2,3,4,5,6,7,8,9,10,11,12]

path='data/paper/paper_fig_2/'
file_pattern = path+"Phase_Diagram_Homophily_numerical_c_4_delta_-12_maxiter_100000_ord_init_conf_0.51_G_*.txt"
file_list = sorted(glob.glob(file_pattern), 
                  key=lambda x: int(x.split('_')[-1].split('.')[0]))  # Sort by G number

# Read each file and store the columns
for file_path in file_list:
    data = np.loadtxt(file_path)
    all_xs_c4.append(data[:, 0])  # First column
    all_ys_c4.append(data[:, 1])  # Second column

all_xs_c4 = np.array(all_xs_c4)
all_ys_c4 = np.array(all_ys_c4)

all_xs_ord_c4 = []
all_ys_ord_c4 = []

file_pattern = path+"Phase_Diagram_Homophily_numerical_c_4_delta_-12_maxiter_100000_ord_init_conf_1.0_G_*.txt"
file_list = sorted(glob.glob(file_pattern), 
                  key=lambda x: int(x.split('_')[-1].split('.')[0]))  # Sort by G number

# Read each file and store the columns
for file_path in file_list:
    data = np.loadtxt(file_path)
    all_xs_ord_c4.append(data[:, 0])  # First column
    all_ys_ord_c4.append(data[:, 1])  # Second column

all_xs_ord_c4 = np.array(all_xs_ord_c4)
all_ys_ord_c4 = np.array(all_ys_ord_c4)



norm_gap_c_4 = [gap_comp_int(all_ys_ord_c4[i], all_ys_c4[i]) if i != 0 else 0 for i in range(len(all_ys_ord_c4))]



substract_gap_c_4=[sum(np.array(all_ys_c4[i])-np.array(all_ys_ord_c4[i])) for i in range(len(all_ys_c4))]

max_gap_c_4=[max(np.array(all_ys_c4[i])-np.array(all_ys_ord_c4[i])) for i in range(len(all_ys_c4))]


integral_c_4=[calculate_area(all_xs_c4[i],all_ys_c4[i]) for i in range(12)]
integral_c_4_ord=[calculate_area(all_xs_ord_c4[i],all_ys_ord_c4[i]) for i in range(12)]

integral_gap_c_4=(np.array(integral_c_4_ord)-np.array(integral_c_4))/np.array(integral_c_4_ord)



def power_law(g, C, alpha):
    return C * (g ** (-alpha))

g_c_4=x_data[2:]
A_c_4=integral_c_4[2:]
# Fit the data
params, covariance = curve_fit(power_law, g_c_4, A_c_4, p0=[1, 1])
C_c_4, alpha_c_4 = params

C_c_4_err = np.sqrt(covariance[0, 0])
alpha_c_4_err = np.sqrt(covariance[1, 1])

print("Power law (c=4):\n")

print(f"Random Fit: A(g) = ({C_c_4:.3f} ± {C_c_4_err:.3f}) * g^(-({alpha_c_4:.3f} ± {alpha_c_4_err:.3f}))")



# Fit the data
params, covariance = curve_fit(power_law, g_c_4, np.array(integral_c_4_ord[2:]), p0=[1, 1])
C_ord_c_4, alpha_ord_c_4 = params

C_ord_c_4_err = np.sqrt(covariance[0, 0])
alpha_ord_c_4_err = np.sqrt(covariance[1, 1])
print(f"Ordered Fit: A(g) = ({C_ord_c_4:.3f} ± {C_ord_c_4_err:.3f}) * g^(-({alpha_ord_c_4:.3f} ± {alpha_ord_c_4_err:.3f}))")


all_xs_c5 = []
all_ys_c5 = []

x_data=[1,2,3,4,5,6,7,8,9,10,11,12]


file_pattern = path+"Phase_Diagram_Homophily_numerical_c_5_delta_-12_maxiter_100000_ord_init_conf_0.51_G_*.txt"
file_list = sorted(glob.glob(file_pattern), 
                  key=lambda x: int(x.split('_')[-1].split('.')[0]))  # Sort by G number

additional_files = [
    path+"Phase_Diagram_Homophily_numerical_c_5_delta_-12_maxiter_100000_rand_init_seed_1_G_11.txt",
    path+"Phase_Diagram_Homophily_numerical_c_5_delta_-12_maxiter_100000_rand_init_seed_1_G_12.txt"
]

# Combine both lists
all_files = file_list + additional_files

# Read each file and store the columns
for file_path in all_files:
    data = np.loadtxt(file_path)
    all_xs_c5.append(data[:, 0])  # First column
    all_ys_c5.append(data[:, 1])  # Second column
    
all_xs_c5 = np.array(all_xs_c5)
all_ys_c5 = np.array(all_ys_c5)

all_xs_ord_c5 = []
all_ys_ord_c5 = []

file_pattern = path+"Phase_Diagram_Homophily_numerical_c_5_delta_-12_maxiter_100000_ord_init_conf_1.0_G_*.txt"
file_list = sorted(glob.glob(file_pattern), 
                  key=lambda x: int(x.split('_')[-1].split('.')[0]))  # Sort by G number

# Read each file and store the columns
for file_path in file_list:
    data = np.loadtxt(file_path)
    all_xs_ord_c5.append(data[:, 0])  # First column
    all_ys_ord_c5.append(data[:, 1])  # Second column

all_xs_ord_c5 = np.array(all_xs_ord_c5)
all_ys_ord_c5 = np.array(all_ys_ord_c5)



norm_gap_c_5 = [gap_comp_int(all_ys_ord_c5[i], all_ys_c5[i]) if i != 0 else 0 for i in range(len(all_ys_ord_c5))]



substract_gap_c_5=[sum(np.array(all_ys_c5[i])-np.array(all_ys_ord_c5[i])) for i in range(len(all_ys_c5))]

max_gap_c_5=[max(np.array(all_ys_c5[i])-np.array(all_ys_ord_c5[i])) for i in range(len(all_ys_c5))]


integral_c_5=[calculate_area(all_xs_c5[i],all_ys_c5[i]) for i in range(12)]
integral_c_5_ord=[calculate_area(all_xs_ord_c5[i],all_ys_ord_c5[i]) for i in range(12)]

integral_gap_c_5=(np.array(integral_c_5_ord)-np.array(integral_c_5))/np.array(integral_c_5_ord)


g_c_5=x_data[2:]
A_c_5=integral_c_5[2:]
# Fit the data
params, covariance = curve_fit(power_law, g_c_5, A_c_5, p0=[1, 1])
C_c_5, alpha_c_5 = params

C_c_5_err = np.sqrt(covariance[0, 0])
alpha_c_5_err = np.sqrt(covariance[1, 1])
print("Power law (c=5):\n")

print(f"Random Fit: A(g) = ({C_c_5:.3f} ± {C_c_5_err:.3f}) * g^(-({alpha_c_5:.3f} ± {alpha_c_5_err:.3f}))")


# Fit the data
params, covariance = curve_fit(power_law, g_c_5, np.array(integral_c_5_ord[2:]), p0=[1, 1])
C_ord_c_5, alpha_ord_c_5 = params

C_ord_c_5_err = np.sqrt(covariance[0, 0])
alpha_ord_c_5_err = np.sqrt(covariance[1, 1])
print(f"Ordered Fit: A(g) = ({C_ord_c_5:.3f} ± {C_ord_c_5_err:.3f}) * g^(-({alpha_ord_c_5:.3f} ± {alpha_ord_c_5_err:.3f}))")



#FIG. 2

plt.figure(figsize=(12, 8), dpi=600) 

main_ax = plt.gca()  

# Main plot formatting (log-log scale)
main_ax.set_xscale('log')
main_ax.set_yscale('log')
main_ax.set_xlabel('g', fontsize=main_xlabel_fs)
main_ax.set_ylabel('A(EC)', fontsize=main_ylabel_fs)
main_ax.tick_params(axis='both', which='major', labelsize=main_xlabel_fs)
main_ax.tick_params(axis='both', which='minor', labelsize=main_ticks_fs)

# Plot power law fits (main plot)
g_fit = np.linspace(3, 12, 100)
main_ax.scatter(g_c_4, integral_c_4[2:], marker='s', s=160, 
                color='darkred', facecolors='darkred', 
                edgecolors='black', zorder=2, label='c=4')
main_ax.plot(g_fit, power_law(g_fit, C_c_4, alpha_c_4), 
             color='red', linestyle='--', zorder=1, linewidth=3.5)

main_ax.scatter(g_c_5, integral_c_5[2:], s=160, 
                color='darkred', facecolors='darkred', 
                edgecolors='black', zorder=2, label='c=5')
main_ax.plot(g_fit, power_law(g_fit, C_c_5, alpha_c_5), 
             color='red', linestyle='--', zorder=1, linewidth=3.5)

# Format log scale ticks
for axis in [main_ax.xaxis, main_ax.yaxis]:
    axis.set_major_formatter(ScalarFormatter(useMathText=True))
    axis.set_minor_formatter(ScalarFormatter())
main_ax.set_xticks([3, 6, 9, 12])
main_ax.set_yticks([0.1, 0.2, 0.3])
main_ax.yaxis.set_major_locator(MaxNLocator(nbins=3)) 

main_ax.set_ylim(0,0.5)

main_ax.minorticks_off()
main_ax.legend(fontsize=main_legend_fs, edgecolor='black', loc='lower left')


# INSET (Gap values) - in bottom-left corner
inset_ax = plt.axes([0.43, 0.48, 0.46, 0.38])  # Position in bottom-left

# Set y-range without extra space
y_max = np.max(np.concatenate([integral_gap_c_4[1:], integral_gap_c_5[1:]]))

# Plot gap data
inset_ax.scatter(x_data[1:], integral_gap_c_4[1:], 
                 color='black', marker='s',
                 facecolors='gray', edgecolors='black',
                 linewidths=1.5, s=70, zorder=2, label='c=4')
inset_ax.plot(x_data[1:], integral_gap_c_4[1:], 
              linewidth=2.2, color='black', zorder=1)

inset_ax.scatter(x_data[1:], integral_gap_c_5[1:], 
                 color='black', marker='o',
                 facecolors='gray', edgecolors='black',
                 linewidths=1.5, s=70, zorder=2, label='c=5')
inset_ax.plot(x_data[1:], integral_gap_c_5[1:], 
              linewidth=2.2, color='black', zorder=1)

# Inset formatting
inset_ax.set_xlabel('g', fontsize=inset_xlabel_fs)
inset_ax.set_ylabel('GS', fontsize=inset_ylabel_fs)
inset_ax.set_xticks([3, 6, 9, 12])

inset_ax.tick_params(axis='both', labelsize=inset_ticks_fs)
inset_ax.legend(fontsize=inset_legend_fs, 
                edgecolor='black', markerscale=1.3, loc='lower right')


# Save and show
plt.savefig('paper_figs/fig_2_gap_c_4_5_powerl_fits_rand.jpg', dpi=800)


#LOADING DATA FIG. 3


file_path='data/paper/paper_fig_3/main/'


file = file_path+'Stability_Para_Homophily_RGfc_c_4_G_3.txt'



# Read the data from the file
data = np.loadtxt(file)

# Separate the data into two arrays
xg3 = data[:, 0]  # First column
yg3 = data[:, 1]  # Second column

# Load data for RANDOM initialization
data_rand = np.loadtxt(file_path+'MC_homophily_Erdos_Renyi_transitions_rand_init_g_3_c_3.txt', skiprows=1)
T_rand = data_rand[:, 0]          # First column: T (beta)
alpha_min_rand = data_rand[:, 1]  # Second column: alpha_min
alpha_max_rand = data_rand[:, 2]  # Third column: alpha_max
alpha_avg_rand = (alpha_min_rand + alpha_max_rand) / 2  # Average alpha

# Load data for ORDERED initialization
data_ord = np.loadtxt(file_path+'MC_homophily_Erdos_Renyi_transitions_ord_init_g_3_c_3.txt', skiprows=1)
T_ord = data_ord[:, 0]            # First column: T (beta)
alpha_min_ord = data_ord[:, 1]    # Second column: alpha_min
alpha_max_ord = data_ord[:, 2]    # Third column: alpha_max
alpha_avg_ord = (alpha_min_ord + alpha_max_ord) / 2  # Average alpha

# Calculate asymmetric errors (distance from mean to min/max)
y_err_lower_rand = alpha_avg_rand - alpha_min_rand
y_err_upper_rand = alpha_max_rand - alpha_avg_rand
y_error_rand = [y_err_lower_rand, y_err_upper_rand]

y_err_lower_ord = alpha_avg_ord - alpha_min_ord
y_err_upper_ord = alpha_max_ord - alpha_avg_ord
y_error_ord = [y_err_lower_ord, y_err_upper_ord]


file = file_path+'BP_ER_popdyn_transitions_c_3_delta_1e-6_maxiter_1000_init_p_0.51_G_3_alphamax_1.0_dalpha_0.01_deltamag_0.01_popsize_100000_smoothsize_100_measuresize_1000_seed_1.txt'



data = np.loadtxt(file)

x_erenyi_c_3_g_3 = data[:, 0]  # First column
y_erenyi_c_3_g_3 = data[:, 1]  # Second column

file = file_path+'BP_ER_popdyn_transitions_c_3_delta_1e-6_maxiter_1000_init_p_1.0_G_3_alphamax_1.0_dalpha_0.01_deltamag_0.01_popsize_100000_smoothsize_100_measuresize_1000_seed_1.txt'



data = np.loadtxt(file)

x_erenyi_c_3_g_3_ord = data[:, 0]  # First column
y_erenyi_c_3_g_3_ord = data[:, 1]  # Second column



file_path='data/paper/paper_fig_3/inset/'
file = file_path+'BP_ER_popdyn_eq_c_3_T_0.200_delta_1e-6_maxiter_1000_init_p_0.51_G_3_deltamag_0.01_popsize_100000_smoothsize_100_measuresize_1000_seed_1.txt'



data = np.loadtxt(file)

f_x = data[:, 1]  # First column
f_y = data[:, 4]  # Second column

file = file_path+'BP_ER_popdyn_eq_c_3_T_0.200_delta_1e-6_maxiter_1000_init_p_1.0_G_3_deltamag_0.01_popsize_100000_smoothsize_100_measuresize_1000_seed_1.txt'



data = np.loadtxt(file)

# Separate the data into two arrays
f_x_ord= data[:, 1]  # First column
f_y_ord = data[:, 4]  # Second column


#FIG. 3



# Create main figure
plt.figure(figsize=(12, 8), dpi=600)

# MAIN PLOT (first plot)
# Plot the main lines
main_ax = plt.gca()  # Get current axes
main_ax.plot(x_erenyi_c_3_g_3, y_erenyi_c_3_g_3, label='BP-rand', color='red', linewidth=3)
main_ax.plot(x_erenyi_c_3_g_3_ord, y_erenyi_c_3_g_3_ord, label='BP-ord', color='blue', linewidth=3)

# Add shaded error regions
main_ax.fill_between(x_erenyi_c_3_g_3, 
                    y_erenyi_c_3_g_3 - 0.01, 
                    y_erenyi_c_3_g_3 + 0.01, 
                    color='red', alpha=0.2)

main_ax.fill_between(x_erenyi_c_3_g_3_ord, 
                    y_erenyi_c_3_g_3_ord - 0.01, 
                    y_erenyi_c_3_g_3_ord + 0.01, 
                    color='blue', alpha=0.2)

# Prediction line
main_ax.plot(xg3, yg3, label='Prediction', linestyle='--', linewidth=3, color='black')

# MC points with error bars
main_ax.errorbar(
    T_rand, alpha_avg_rand, 
    yerr=y_error_rand, fmt='s', markersize=10,
    markerfacecolor='red', color='black', label='MC-rand',
    capsize=5, elinewidth=2
)
main_ax.errorbar(
    T_ord, alpha_avg_ord, 
    yerr=y_error_ord, 
    fmt='^', markersize=10,
    markerfacecolor='blue', color='black', label='MC-ord',
    capsize=5, elinewidth=2
)

main_ax.legend(fontsize=main_legend_fs, loc='lower right')
main_ax.set_ylabel(r'$\alpha_c$', fontsize=main_ylabel_fs)
main_ax.set_xlabel('T', fontsize=main_xlabel_fs)
main_ax.tick_params(axis='both', labelsize=main_ticks_fs)

# CREATE INSET (second plot)
inset_ax = inset_axes(main_ax, 
                     width="60%", 
                     height="60%",
                     loc='upper left',
                     bbox_to_anchor=(0.10, 0.3, 0.7, 0.7),
                     bbox_transform=main_ax.transAxes)

# Plot inset content (second plot)
inset_ax.scatter(f_x, f_y, 
                color='red', label='BP-rand', marker='s', 
                edgecolor='red', facecolor='none', linewidths=1, s=25)
inset_ax.scatter(f_x_ord, f_y_ord, 
                color='blue', label='BP-ord', marker='^', 
                edgecolor='blue', facecolor='none', linewidths=1, s=25)

inset_ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
inset_ax.legend(fontsize=inset_legend_fs, markerscale=3)  # Slightly reduced font size
inset_ax.set_ylabel('f', fontsize=inset_ylabel_fs)  # Slightly reduced font size
inset_ax.set_xlabel(r'$\alpha$', fontsize=inset_xlabel_fs)  # Slightly reduced font size
inset_ax.tick_params(axis='both', labelsize=inset_ticks_fs)  # Slightly reduced tick size
inset_ax.text(0.82, -6, 'T = 0.2', rotation=0, fontsize=26, color='black')  # Slightly reduced font size


plt.savefig('paper_figs/fig_3_er_pd_free_e_T_0.2.jpg', dpi=800)

