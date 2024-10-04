include("./Optimizing_t_Estimator_commented.jl")

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

#GRAPH

graphname = ARGS[9]
graphid = ARGS[10]

init_conf = ARGS[11]

g = create_graph(n, graphname)

edge_n = fill_edge_n(g)


fpl="plinks_status"
fcor="correlations"
fme="Mean_Energy"
fmago="Mag_ord"
fmagd="Mag_deord"


all(edge_n, t, ts, g, G, mcsteps, n, k, N, αv, βv, fpl, fcor, fme, fmago, fmagd, graphid, init_conf)