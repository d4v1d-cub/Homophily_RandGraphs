import numpy as np
from glob import glob

def detect_transition(filename, delta_mag, dalpha):
    try:
        fin = open(filename, "r")
        while True:
            j = fin.readline()
            if not j:
                return -1, -1  # End of file reached
            if not j.startswith("#"):
                line = j.split()
                temp = float(line[0])
                alpha = float(line[1])
                mag = float(line[3])
                if mag > delta_mag:
                    return temp, alpha - dalpha / 2
    except (IOError, OSError) as error:
        print("An error occurred: ", error)
        return -1, -1
    

def read_all(path, string1, string2, string3, delta_mag, dalpha):
    trans_alphas = []
    temps = []
    filelist = glob(f'{path}/AllData/{string1}**{string2}**{string3}')
    for filename in filelist:
        temp, alpha = detect_transition(filename, delta_mag, dalpha)
        if temp > 0:
            temps.append(temp)
            trans_alphas.append(alpha)
    sorted_temps_alphas = sorted(zip(temps, trans_alphas), key=lambda x: x[0])
    return sorted_temps_alphas


def print_transition(path, fileout, sorted_temps_alphas):
    try:
        with open(f'{path}/{fileout}', "w") as fout:
            fout.write("# T alpha\n")
            for temp, alpha in sorted_temps_alphas:
                fout.write(f"{temp}\t{alpha}\n")
    except (IOError, OSError) as error:
        print("An error occurred while writing to file: ", error)



cm = "3"
delta = "1e-6"
maxiter = "1000"
init_p = "1.0"
G = "3"
alphamax = "1.0"
dalpha = 0.01
popsize = "100000"
smoothsize = "100"
measuresize = "1000"
seed = "1"

delta_mag = 0.01


pathtofiles = f'/media/david/Data/UH/Grupo_de_investigacion/Homophily/BP/RGER/Results/'
string1 = f'BP_ER_popdyn_eq_c_{cm}_T_'
string2 = f'_delta_{delta}_maxiter_{maxiter}_init_p_{init_p}_G_{G}_alphamin_'
string3 = f'_alphamax_{alphamax}_dalpha_{dalpha}_deltamag_{delta_mag}_popsize_{popsize}_smoothsize_{smoothsize}_measuresize_{measuresize}_seed_{seed}.txt'

fileout = f'BP_ER_popdyn_transitions_c_{cm}_delta_{delta}_maxiter_{maxiter}_init_p_{init_p}_G_{G}_alphamax_{alphamax}_dalpha_{dalpha}_deltamag_{delta_mag}_popsize_{popsize}_smoothsize_{smoothsize}_measuresize_{measuresize}_seed_{seed}.txt'

sorted_temps_alphas = read_all(pathtofiles, string1, string2, string3, delta_mag, dalpha)
if sorted_temps_alphas:
    print_transition(pathtofiles, fileout, sorted_temps_alphas)