using LightGraphs, LinearAlgebra 
import Random


# This function normalizes the messages after they are computed using the sum-product rule
function normalize_norm_1(vec::Array{Float64, 1})
    norm = sum(vec)
    for i in 1:length(vec)
        vec[i] /= norm
    end
    return vec
end


# This function initalizes the graph as an instance of a LightGraphs object.
# The information about the graph is read from the file with name "filename". The graph size is "n"
# The file has the following format. For each node you have
# number_of_neighbors  ----->   Integer
# node_index     neighbor_index        Jij      Jji
function Import_graph_from_file(filename::String, n::Int64)
    graph = Graph(n)
    r = open(filename, "r")
    node = 1
    while !eof(r)
        numneigh = parse(Int64, readline(r))
        for i in 1:numneigh
            line = split(readline(r))
            neigh = parse(Int64, line[2]) + 1     # As julia takes indexes starting from 1, I add 1 to the index of the neighbor
            add_edge!(graph, node, neigh)
        end
        node += 1
    end
    close(r)
    return graph
end



# This function initializes all messages. The messages will be stores in a bi-dimensional 
# Array{Array{Float64, 1}, 1}. The first index goes twice over all possible edges in the graph.
# Twice because there are two directions associated to each edge. The message m_i->j and the 
# message m_j->i. 
# The second index goes over all the messages associated to one edge that go in the same direction.
# There are 2^{G} different messages m_i->j, because the vector sitting in the site j can have
# 2^{G} different values. 
function set_messages(graph::Graph, seed::Int64, G::Int64)
    Random.seed!(seed)
    dict_messages = Dict{String, Int64}()   # This dictionary says what is the place in the arrays
                                            # "messages" and "messages_new" of the message m_i_j
    messages = Array{Array{Float64, 1}, 1}()
    messages_new = Array{Array{Float64, 1}, 1}()     # This is meant to temporarily store the updates 
                                                     # of the messages during a sweep of the sum-product 
                                                     # update rule 
    counter = 1
    for edge in edges(graph)
        mess = Array{Float64, 1}()
        mij = "m_" * string(src(edge)) * "_" * string(dst(edge))  # The keys of the dictionary are strings
        dict_messages[mij] = counter
        for comb in 1:2^G
            rnd = rand()              # At the beginning, the messages are random
            push!(mess, rnd)
        end
        mess = normalize_norm_1(mess)    # I normalize
        push!(messages, deepcopy(mess))     # deepcopy is needed to instantiate two different arrays
        push!(messages_new, deepcopy(mess))

        counter += 1

        mess = Array{Float64, 1}()
        mji = "m_" * string(dst(edge)) * "_" * string(src(edge))
        dict_messages[mji] = counter
        for comb in 1:2^G
            rnd = rand()
            push!(mess, rnd)
        end
        mess = normalize_norm_1(mess)
        push!(messages, deepcopy(mess))
        push!(messages_new, deepcopy(mess))     
    end
    return messages, messages_new, dict_messages
end


function prod_others(graph::Graph, messages::Array{Array{Float64, 1}, 1}, 
                     dict_messages::Dict{String, Int64}, nodo::Int64, neigh::Int64, comb::Int64)
    prod = 1.0
    others = filter(x -> x != neigh , neighbors(graph, nodo))
    for k in others
        place = dict_messages["m_" * string(k) * "_" * string(nodo)]
        prod *= messages[place][comb]    # Product over the incoming messages m_k->j(s_j)
    end
    return prod
end

# This is the function that computes the new messages
# spins_combs[comb] gives a vector with elements 1 and -1, that corresponds to the integer "comb"
function sum_prod_rule(graph::Graph, temp::Float64, G::Int64, alpha::Float64, 
                       spins_combs::Array{Array{Int8, 1}, 1}, messages::Array{Array{Float64, 1}, 1}, 
                       messages_new::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64})
    for edge in edges(graph)
        for (i, j) in [(src(edge), dst(edge)), (dst(edge), src(edge))]
            place_ij = dict_messages["m_" * string(j) * "_" * string(i)]
            for comb_i in 1:2^G
                cumulant = 0.0
                for comb_j in 1:2^G
                    prod = prod_others(graph, messages, dict_messages, j, i, comb_j)
                    d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])   # This is s_i . s_j
                    prod *= exp(((2 * alpha - 1) * d_prod + abs(d_prod)) / temp / 2 / G)
                    # The pair interaction has the shape (2 α - 1) s_i . s_j / 2 G  +  |s_i . s_j| / 2 G
                    # Notice that if α=0, the sum of this two terms is zero for s_i . s_j > 0, and greater 
                    # than zero if s_i . s_j < 0. Therefore, the heterophylic relations are favored
                    # If α=1, the opposite happens and homophylic relations are favored
                    cumulant += prod
                end
                messages_new[place_ij][comb_i] = cumulant
            end
            messages_new[place_ij] = normalize_norm_1(messages_new[place_ij])
        end
    end
    return messages_new
end


# This function computes a (defined by me) distance between the old messages "messages" and the 
# updated messages "messages_new". distance(messages, messages_new) = 
# max(|messages - messages_new| / messages)
# I now realize that I could just code that directly in julia array-like operations. This is 
# old coding.
function comp_error(graph::Graph, G::Int64, messages::Array{Array{Float64, 1}, 1}, 
                    messages_new::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64})
    max_delta = -1
    for edge in edges(graph)
        for (i, j) in [(src(edge), dst(edge)), (dst(edge), src(edge))]
            for comb_i in 1:2^G
                place = dict_messages["m_" * string(j) * "_" * string(i)]
                delta = abs(messages_new[place][comb_i] - messages[place][comb_i]) /
                            messages[place][comb_i]
                messages[place][comb_i] = messages_new[place][comb_i]
                if delta > max_delta
                    max_delta = delta
                end
            end
        end
    end
    return max_delta, messages
end


# This is self-explanatory
function update_messages(graph::Graph, temp::Float64, G::Int64, alpha::Float64, 
                         spins_combs::Array{Array{Int8, 1}, 1}, messages::Array{Array{Float64, 1}, 1}, 
                         messages_new::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64})
    
    messages_new = sum_prod_rule(graph, temp, G, alpha, spins_combs, messages, messages_new, 
                                 dict_messages)
    error, messages = comp_error(graph, G, messages, messages_new, dict_messages)
    return error
end


# Computes the average dot product <s_i . s_j>
function ener(graph::Graph, G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
              messages::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64}, 
              alpha::Float64, temp::Float64)
    cumulant = 0.0
    for edge in edges(graph)
        cumul_edge = 0
        i = src(edge)
        j = dst(edge)
        cumul_prob = 0.0
        for comb_i in 1:2^G
            for comb_j in 1:2^G
                prodi = prod_others(graph, messages, dict_messages, i, j, comb_i)
                prodj = prod_others(graph, messages, dict_messages, j, i, comb_j)
                d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])
                hamilt = ((2 * alpha - 1) * d_prod + abs(d_prod)) / 2 / G
                exp_hamilt = exp(hamilt / temp)
                cumul_edge += prodi * prodj * hamilt * exp_hamilt
                cumul_prob += prodi * prodj * exp_hamilt
            end
        end
        cumulant += cumul_edge / cumul_prob
    end
    return cumulant / ne(graph)
end


# Computes the average dot product <s_i . s_j>
function av_pair(graph::Graph, G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
                 messages::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64},
                 alpha::Float64, temp::Float64)
    cumulant = 0.0
    for edge in edges(graph)
        cumul_edge = 0
        i = src(edge)
        j = dst(edge)
        cumul_prob = 0.0
        for comb_i in 1:2^G
            for comb_j in 1:2^G
                prodi = prod_others(graph, messages, dict_messages, i, j, comb_i)
                prodj = prod_others(graph, messages, dict_messages, j, i, comb_j)
                d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])
                exp_hamilt = exp(((2 * alpha - 1) * d_prod + abs(d_prod)) / temp / 2 / G)
                cumul_edge += prodi * prodj * d_prod / G * exp_hamilt
                cumul_prob += prodi * prodj * exp_hamilt
            end
        end
        cumulant += cumul_edge / cumul_prob
    end
    return cumulant / ne(graph)
end



function av_abs_pair(graph::Graph, G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
                 messages::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64},
                 alpha::Float64, temp::Float64)
    cumulant = 0.0
    for edge in edges(graph)
        cumul_edge = 0
        i = src(edge)
        j = dst(edge)
        cumul_prob = 0.0
        for comb_i in 1:2^G
            for comb_j in 1:2^G
                prodi = prod_others(graph, messages, dict_messages, i, j, comb_i)
                prodj = prod_others(graph, messages, dict_messages, j, i, comb_j)
                d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])
                exp_hamilt = exp(((2 * alpha - 1) * d_prod + abs(d_prod)) / temp / 2 / G)
                cumul_edge += prodi * prodj * abs(d_prod) / G * exp_hamilt
                cumul_prob += prodi * prodj * exp_hamilt
            end
        end
        cumulant += cumul_edge / cumul_prob
    end
    return cumulant / ne(graph)
end


# Computes the average probability of having a positive interaction s_i . s_j > 0
function pos_links(graph::Graph, G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
                   messages::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64},
                   alpha::Float64, temp::Float64)
    cumulant = 0.0
    for edge in edges(graph)
        cumul_edge = 0
        i = src(edge)
        j = dst(edge)
        cumul_prob = 0.0
        for comb_i in 1:2^G
            for comb_j in 1:2^G
                prodi = prod_others(graph, messages, dict_messages, i, j, comb_i)
                prodj = prod_others(graph, messages, dict_messages, j, i, comb_j)
                d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])
                exp_hamilt = exp(((2 * alpha - 1) * d_prod + abs(d_prod)) / temp / 2 / G)
                if d_prod > 0
                    # I only sum if s_i . s_j > 0
                    cumul_edge += prodi * prodj * exp_hamilt
                end
                cumul_prob += prodi * prodj * exp_hamilt
            end
        end
        cumulant += cumul_edge / cumul_prob
    end
    return cumulant / ne(graph)
end



function neg_links(graph::Graph, G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
                   messages::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64},
                   alpha::Float64, temp::Float64)
    cumulant = 0.0
    for edge in edges(graph)
        cumul_edge = 0
        i = src(edge)
        j = dst(edge)
        cumul_prob = 0.0
        for comb_i in 1:2^G
            for comb_j in 1:2^G
                prodi = prod_others(graph, messages, dict_messages, i, j, comb_i)
                prodj = prod_others(graph, messages, dict_messages, j, i, comb_j)
                d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])
                exp_hamilt = exp(((2 * alpha - 1) * d_prod + abs(d_prod)) / temp / 2 / G)
                if d_prod < 0
                    # I only sum if s_i . s_j > 0
                    cumul_edge += prodi * prodj * exp_hamilt
                end
                cumul_prob += prodi * prodj * exp_hamilt
            end
        end
        cumulant += cumul_edge / cumul_prob
    end
    return cumulant / ne(graph)
end


function zero_links(graph::Graph, G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
                    messages::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64},
                    alpha::Float64, temp::Float64)
    cumulant = 0.0
    for edge in edges(graph)
        cumul_edge = 0
        i = src(edge)
        j = dst(edge)
        cumul_prob = 0.0
        for comb_i in 1:2^G
            for comb_j in 1:2^G
                prodi = prod_others(graph, messages, dict_messages, i, j, comb_i)
                prodj = prod_others(graph, messages, dict_messages, j, i, comb_j)
                d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])
                exp_hamilt = exp(((2 * alpha - 1) * d_prod + abs(d_prod)) / temp / 2 / G)
                if d_prod == 0
                    # I only sum if s_i . s_j > 0
                    cumul_edge += prodi * prodj * exp_hamilt
                end
                cumul_prob += prodi * prodj * exp_hamilt
            end
        end
        cumulant += cumul_edge / cumul_prob
    end
    return cumulant / ne(graph)
end


function global_mag(graph::Graph, G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
    messages::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64})
    cumul_prob_i = zeros(2^G)
    prob_i = ones(2^G)
    for node in vertices(graph)
        for comb_i in 1:2^G
            prod = 1.0
            neighs = neighbors(graph, node)
            for k in neighs
                place = dict_messages["m_" * string(k) * "_" * string(node)]
                prod *= messages[place][comb_i]    # Product over the incoming messages m_k->i(s_i)
            end
            prob_i[comb_i] = prod
        end
        prob_i /= sum(prob_i)

        cumul_prob_i += prob_i
    end

    cumul_prob_i /= nv(graph)
    index_max = findmax(cumul_prob_i)[2]
    cumul = 0.0
    for comb_i in 1:2^G
        norm_d_prod = dot(spins_combs[comb_i], spins_combs[index_max]) / G
        cumul += norm_d_prod * cumul_prob_i[comb_i]
    end 
    return cumul
end


function local_mag(graph::Graph, G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
    messages::Array{Array{Float64, 1}, 1}, dict_messages::Dict{String, Int64})
    cumulant = 0.0
    prob_i = ones(2^G)
    for node in vertices(graph)
        for comb_i in 1:2^G
            prod = 1.0
            neighs = neighbors(graph, node)
            for k in neighs
                place = dict_messages["m_" * string(k) * "_" * string(node)]
                prod *= messages[place][comb_i]    # Product over the incoming messages m_k->i(s_i)
            end
            prob_i[comb_i] = prod
        end
        index_max = findmax(prob_i)[2]
        prob_i /= sum(prob_i)

        cumul_node = 0.0
        for comb_i in 1:2^G
            norm_d_prod = dot(spins_combs[comb_i], spins_combs[index_max]) / G
            cumul_node += norm_d_prod * prob_i[comb_i]
        end
        cumulant += cumul_node
    end
    return cumulant / nv(graph)
end


# This funtion converts an integer into an array of elements 1 and -1.
function int_to_spins(comb::Int64, G::Int64)
    str = bitstring(comb)   # The integer is encoden into a string of 0 and 1 (bits)
    spins = Array{Int8, 1}()
    for i in 0:G-1
        push!(spins, 1 - 2 * parse(Int64, str[end-i]))   # The array "spins" is filled with 1 and -1
                                                         # When the bit is 0, the spin is 1
                                                         # When the bit is 1, the spin is -1
    end
    return spins
end


# This function creates an Array{Array{Int8, 1}, 1} called "spins_combs"
# spins_combs[comb] gives the corresponding array of 1 and -1
function get_spins_combs(G::Int64)
    spins_combs = Array{Array{Int8, 1}, 1}()
    for comb in 1:2^G
        push!(spins_combs, int_to_spins(comb - 1, G))
    end
    return spins_combs
end


# This function runs a single BP convergence
function belief_propagation(graph::Graph, temp::Float64, max_iter::Int64, sampling::Int64,
    G::Int64, alpha::Float64, spins_combs::Array{Array{Int8, 1}, 1}, 
    messages::Array{Array{Float64, 1}, 1}, messages_new::Array{Array{Float64, 1}, 1}, 
    dict_messages::Dict{String, Int64}, delta::Float64)
    
    iter = 0
    max_delta = delta + 1
    println("T=", temp, "\t", "alpha=", alpha, "\t", "iter=", iter, "\t", "error=", max_delta)
    k = 1
    while (iter < max_iter && max_delta > delta)
        max_delta = update_messages(graph, temp, G, alpha, spins_combs, messages, messages_new, 
                                    dict_messages)
        if (iter == k * sampling)
            println("T=", temp, "\t", "alpha=", alpha, "\t", "iter=", iter, "\t", "error=", max_delta)
            k += 1
        end
        iter += 1
    end
    println("T=", temp, "\t", "alpha=", alpha, "\t", "iter=", iter, "\t", "error=", max_delta)
    return max_delta < delta
end


# This function runs several BP convergences for different parameters. In this case, the parameters 
# are T and α
function run_BP_several_params(filename, seed, string1, n, max_iter, sampling, G, alphas_list, delta)
    spins_combs = get_spins_combs(G)
    graph = Import_graph_from_file(filename, n)
    for temp in temps
        fileeq = string("BP_", string1,"_eq_N_", string(n), "_c_", string(c), "_T_", string(temp), 
                        "_delta_", string(exponent), "_sampling_", string(sampling), "_maxiter_", 
                        string(max_iter), "_seed_", string(seed), "_G_", string(G), ".txt")
        w = open(fileeq, "w") 
        for alpha in alphas_list         
            messages, messages_new, dict_messages = set_messages(graph, seed, G)
            conv = belief_propagation(graph, temp, max_iter, sampling, G, alpha, spins_combs, messages, 
                                      messages_new, dict_messages, delta)
            eq_e = ener(graph, G, spins_combs, messages, dict_messages, alpha, temp)
            eq_av_dot = av_pair(graph, G, spins_combs, messages, dict_messages, alpha, temp)
            eq_abs_dot = av_abs_pair(graph, G, spins_combs, messages, dict_messages, alpha, temp)
            eq_pl = pos_links(graph, G, spins_combs, messages, dict_messages, alpha, temp)
            eq_nl = neg_links(graph, G, spins_combs, messages, dict_messages, alpha, temp)
            eq_zl = zero_links(graph, G, spins_combs, messages, dict_messages, alpha, temp)
            eq_gl_m = global_mag(graph, G, spins_combs, messages, dict_messages)
            eq_loc_m = local_mag(graph, G, spins_combs, messages, dict_messages)
            println(temp, "\t", alpha, "\t", "e=", eq_e, "\t", "avdot=", eq_av_dot, "\t", 
                    "absdot=", eq_abs_dot, "\t", "pos_l=", eq_pl, "\t", "neg_l=", eq_nl, "\t", 
                    "zero_l=", eq_zl, "\t", "global m=", eq_gl_m, "\t", "local m=", eq_loc_m, "\t", 
                    "converged=", conv)
            write(w, string(temp) * "\t" * string(alpha) * "\t" * string(eq_e) * 
                     "\t" * string(eq_av_dot) * "\t" * string(eq_abs_dot) * 
                     "\t" * string(eq_pl) * "\t" * string(eq_nl) * 
                     "\t" * string(eq_zl) * "\t" * string(eq_gl_m) * "\t" * string(eq_loc_m) * 
                     "\t" * string(conv) * "\n")
        end
        close(w)
    end
end


n = 1000   # Number of nodes in the graph
c = 3      # Graph connectivity
idum = -2  # This is the seed of the graph that you will read from a graph
# You can change this part if you want to see what happens with other graphs. The only thing that 
# you need is to respect the format that I use for storing the information inside the files.
filename = string("RanGraph_N_", string(n), "_fixed_c_", string(c), "_simetric_1_model_1_idum1_", string(idum), ".txt")

temps = [0.2, 0.5]
exponent = -15
delta = 10.0 ^ exponent   # When the error is below this threshold, I declare that BP converges
sampling = 10             # BP will print output every "sampling" iterations
max_iter = 1000            # BP will stop after "max_iter" iterations
seed = 1                  # The seed for the initial messages
string1 = "RGfc"          # An identification of the graph you are using. This string is insterted in 
                          # the output strings
G = parse(Int64, ARGS[1])   
# G = 2               
alphas_list = range(0.5, 1, step=0.05)

run_BP_several_params(filename, seed, string1, n, max_iter, sampling, G, alphas_list, delta)