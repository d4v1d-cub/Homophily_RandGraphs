using Erdos

n=parse(Int64, ARGS[1])
k=parse(Int64, ARGS[2])
graphname=ARGS[3]
seed=parse(Int64, ARGS[4])

if graphname=="Barabasi_Albert"
	g=barabasi_albert(n, k, seed=seed)
elseif graphname=="Erdos_Renyi"
	m = n * k ÷ 2
	g=erdos_renyi(n, m, seed=seed)
elseif graphname=="Random_Regular"
	g=random_regular_graph(n, k, seed=seed)
end


function graph_export(graphid,n,k,seed)
    open("Graph_$(graphid)_N_$(n)_k_$(k)_seed_$(seed).txt", "w") do f
        m = n * k ÷ 2
        write(f, "$(n) $(m)\n")
        for e in edges(g)
            write(f, "$(src(e)) $(dst(e))\n")
        end
    end
end
    
graph_export(graphname,n,k,seed)

