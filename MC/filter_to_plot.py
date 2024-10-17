import numpy as np
from glob import glob

def parse_file(filename, nsteps, dt, fileout):
    try:
        fin = open(filename, "r")
        av = np.zeros(nsteps + 1)
        av_sqr = np.zeros(nsteps + 1)
        counter = 0
        nsampl = np.zeros(nsteps + 1)
        while True:
            j = fin.readline()
            if not j:
                break
            nsampl[counter] += 1
            av[counter] += float(j)
            av_sqr[counter] += float(j) ** 2
            counter = (counter + 1) % (nsteps + 1)
        av /= nsampl
        av_sqr /= nsampl
        fin.close()
        fout = open(fileout, "w")
        times = np.arange(0, nsteps * dt + dt / 2, dt)
        
        for i in range(nsteps):
            error = np.sqrt(abs(av_sqr[i] - av[i] ** 2) / nsteps)
            fout.write(f'{times[i]}\t{av[i]}\t{error}\n')
        fout.close()
    except (IOError, OSError) as error:
        print("An error occurred: ", error)


obs = "Mag_deord"
n = 10000
c = 4
g = 4
beta_list = ["2.0", "2.13", "2.5", "3.03", "3.85", "5.26"]
sims = 10
nsteps = 100
dt = 1000
cond_init = 'ord'

pathtofiles = f'/media/david/Data/UH/Grupo_de_investigacion/Homophily/MC/Results/Transition/N_{n}'

for beta in beta_list:
    filelist = glob(f'{pathtofiles}/Full_Random_Regular_{obs}_beta_{beta}_nodes_{n}_neigh_{c}_G_{g}_MCsteps_{nsteps * dt}_sims_{sims}_alpha_**_{cond_init}.dat')
    for filein in filelist:
        filename = filein.split("/")[-1]
        print(f'Reading file: "{filename}"')
        fileout = pathtofiles + f'/Parsed/G_{g}/' + filename[:-4] + "_parsed.txt"
        parse_file(filein, nsteps, dt, fileout)