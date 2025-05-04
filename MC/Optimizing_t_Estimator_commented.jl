using Erdos
using LinearAlgebra
using DelimitedFiles
using JSON


#FUNCTIONS


function opinion_vector(G)  #Gen the opinion vector
    v=rand(0:1,G)
    for i=1:G
        if v[i]==0
            v[i]=-1
        end
    end
    return v
end


function opinion_vector_o(G)  #Gen the opinion vector
    v=ones(G)
    return v
end



function weight_comp(g, Global_State) #Computing the weights
    return reduce(vcat,[[dot(Global_State[i], Global_State[v]) for v in neighbors(g,i)] for i=1:nv(g)])
end


function init(g,G,n, init_conf)
    
    
    Global_State=[]

    if init_conf == "ord"
        for _=1:n
            push!(Global_State, opinion_vector_o(G))
        end
    elseif init_conf == "rand"
        for _=1:n
            push!(Global_State, opinion_vector(G))
        end
    else
        println("init_conf must be 'ord' or 'rand'")
        exit(1)
    end
    
    weights_v = weight_comp(g, Global_State)

    return weights_v, Global_State
end



function stress(g, knode, α, G, weights_v, edge_n)
    friendship = 0.0
    rivarly = 0.0
    for it in eachindex(neighbors(g, knode))
        w=weights_v[edge_n[knode][it]]
        if w > 0
            friendship += w
        elseif w < 0
            rivarly += w
        end
    end
    return -(α / G) * friendship + ((1 - α) / G) * rivarly
end


function stress_2(α, G,weight_new)
    friendship = 0.0
    rivarly = 0.0
    for it in eachindex(weight_new)
        w=weight_new[it]
        if w > 0
            friendship += w
        elseif w < 0
            rivarly += w
        end
    end
    return -(α / G) * friendship + ((1 - α) / G) * rivarly
end


function mean_ener(g, α, G, weights_v, edge_n)
    s = Vector{Float64}(undef, n)
    for i in 1:n
        s[i] = stress(g, i, α, G, weights_v, edge_n)
    end
    return sum(s) / n
end


function mag_o(G,n, Global_State)
    s0=ones(G)
    s = Vector{Float64}(undef, n)
    for i in 1:n
        s[i] = dot(s0,Global_State[i])
    end
    return sum(s)/(G*n)
end

function mag_d(G,n, Global_State)
    return norm(sum(Global_State))/(n*sqrt(G))
end

function update_weights(kn,wn,wv,en)
    for i in eachindex(wn)
        wv[en[kn][i]]=deepcopy(wn[i])
    end
end


function MC_step(g, G, n, α, β, weights_v, edge_n, Global_State)
    
    knode = rand(1:n)
    h1 = stress(g, knode, α, G, weights_v, edge_n)
    r = rand(1:G)
    
    Global_State_knode_r = deepcopy(Global_State[knode][r])
    Global_State[knode][r] *= -1
    
    neighbors_i = neighbors(g, knode)
    weight_new= [deepcopy(weights_v[edge_n[knode][it]])+2*Global_State[knode][r]*Global_State[j][r] for (it,j) in  enumerate(neighbors_i)] 
    
    h2 = stress_2(α, G, weight_new)
    exp_val = exp(-β*(h2-h1))  
    
    if rand() < exp_val
        weights_v = weight_comp(g, Global_State)
    else
        Global_State[knode][r] = Global_State_knode_r
    end
    return weights_v, Global_State
end


function plinks_count(g, weights_v) #Counting positive links
    plc=0
    for w in weights_v
        if w>0
            plc+=1
        end
    end
    return plc/(2*ne(g))
end

function correlation(g,G, weights_v) 
    corr=0
    for w in weights_v
        corr+=w    
    end
    return corr/(G*2*ne(g))
end


function save_data(A,β,α,fn,n,k,N,G,mcsteps,graphid, init_conf)
    io=open("Results/Full_$graphid"*"_$fn"*"_beta_$β"*"_nodes_$n"*"_neigh_$k"*"_G_$G"*"_MCsteps_$mcsteps"*"_sims_$N"*
            "_alpha_$α"*"_"*init_conf*".dat", "a")
    writedlm(io, A)
    close(io)
end


#function save_data_w(A,β,fn)
#    io=open("Results/Full_$fn"*"_$β.dat", "a")
#    for v in A
#        data=JSON.json(v)
#        write(io, data,"\n")
#    end
#    close(io)
#end

#function annotate_data(fn,β,α)
#    io=open("Results/Full_$fn"*"_$β.dat", "a")
#        writedlm(io, α,"\n")
#    close(io)
#end


function main(edge_n, weights_v, Global_State, t, ts, g, G, mcsteps, n, k, N, α, β, fpl, fcor, 
              fme, fmago, fmagd, graphid, init_conf) #Core function

    plinks_status=[]
    me=[]
    mago=[]
    magd=[]
    corr_vector=[]

    for i=1:t
        weights_v, Global_State = MC_step(g,G,n,α,β, weights_v, edge_n, Global_State)
        if i%(ts)==0 || i==1
            #println("t="," ","$i")
            push!(me, mean_ener(g, α,G, weights_v, edge_n))
            push!(plinks_status, plinks_count(g, weights_v))
            push!(corr_vector,correlation(g, G, weights_v))
            push!(mago, mag_o(G,n, Global_State))
            push!(magd, mag_d(G,n, Global_State))
        end
    end
    save_data(plinks_status,β,α,fpl,n,k,N,G,mcsteps,graphid, init_conf)
    save_data(corr_vector,β,α,fcor,n,k,N,G,mcsteps, graphid, init_conf)
    save_data(me,β,α,fme,n,k,N,G,mcsteps, graphid, init_conf)
    save_data(mago,β,α,fmago,n,k,N,G,mcsteps, graphid, init_conf)
    save_data(magd,β,α,fmagd,n,k,N,G,mcsteps, graphid, init_conf)
end

#VARS

function all(edge_n, t, ts, g, G, mcsteps, n, k, N, αv, βv, fpl, fcor, fme, fmago, fmagd, graphid, 
             init_conf)
    for l=1:N
        weights_v, Global_State = init(g,G,n, init_conf)
        main(edge_n, weights_v, Global_State, t, ts,g,G,mcsteps,n,k,N,αv,βv, fpl, fcor, 
             fme, fmago, fmagd, graphid, init_conf)
        println("trajectory no. $l  done")
    end
end


function create_graph(n, graphname)
    g = Graph(n)
    open(graphname, "r") do f
        readline(f) # Skip the first line
        for line in eachline(f)
            src, dst = split(line)
            add_edge!(g, parse(Int, src), parse(Int, dst))
        end
    end
    return g
end


function fill_edge_n(g)
    cumulneighs = 0
    edge_n = []
    for i in 1:nv(g)
        push!(edge_n, (1:length(neighbors(g,i))) .+ cumulneighs)
        cumulneighs += length(neighbors(g,i))
    end
    return edge_n
end