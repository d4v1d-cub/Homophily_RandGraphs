include("./Optimizing_t_Estimator_commented.jl")

# αv=0.8
# βv=1.0

# n=10000
# k=6
# G=9
# mcsteps=1
# N=2

# t=n*G*mcsteps


# graphid = "ER"
# graphstr0 = "Graph_Erdos_Renyi_N_10000_k_6"

# gseed0 = 1
# gstep = 1
# gseedf = 2


αv=parse(Float32, ARGS[1])
βv=parse(Float32, ARGS[2])

n=parse(Int64, ARGS[3])
k=parse(Int64, ARGS[4])
G=parse(Int64, ARGS[5])
mcsteps=parse(Int64, ARGS[6])
tsampling=parse(Int64, ARGS[7])
N=parse(Int64, ARGS[8])

t=n*G*mcsteps
ts=n*G*tsampling

graphid = ARGS[9]
graphstr0 = ARGS[10]

gseed0 = parse(Int64, ARGS[11])
gstep = parse(Int64, ARGS[12])
gseedf = parse(Int64, ARGS[13])

init_conf = ARGS[14]

for seed in gseed0:gstep:gseedf
    graphname = graphstr0 * "_seed_$seed.txt"
    g = create_graph(n, graphname)

    edge_n = fill_edge_n(g)    

    fpl="plinks_status"
    fcor="correlations"
    fme="Mean_Energy"
    fmago="Mag_ord"
    fmagd="Mag_deord"


    all(edge_n, t, ts, g, G, mcsteps, n, k, N, αv, βv, fpl, fcor, fme, fmago, fmagd, graphid, init_conf)

    println("Graph seed $seed  done\n")
end