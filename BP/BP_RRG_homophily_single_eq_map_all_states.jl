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

function set_messages_rand(seed::Int64, G::Int64)
    Random.seed!(seed)

    messages = Array{Float64, 1}()
    messages_new = Array{Float64, 1}()   
    mess = Array{Float64, 1}()
    for _ in 1:2^G
        rnd = rand()              # At the beginning, the messages are random
        push!(mess, rnd)
    end
    mess = normalize_norm_1(mess)    
    messages=deepcopy(mess)
    messages_new=deepcopy(mess)
    return messages, messages_new
end


function set_messages_ord(conf::Int64, G::Int64)
    messages=zeros(2^G)
    messages[conf] = 1
    messages_new=deepcopy(messages)
    return messages, messages_new
end


function sum_prod_rule(temp::Float64, c::Int64, G::Int64, alpha::Float64, 
                       spins_combs::Array{Array{Int8, 1}, 1}, messages::Array{Float64, 1}, 
                       messages_new::Array{Float64, 1})
        for comb_i in 1:2^G
            cumulant = 0.0
            for comb_j in 1:2^G
                d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])   # This is s_i . s_j
                prod = messages[comb_j] * exp(((2 * alpha - 1) * d_prod + abs(d_prod)) / temp / 2 / G) 
                cumulant += prod
            end
            messages_new[comb_i] = cumulant^(c-1)
        end
        messages_new = normalize_norm_1(messages_new)
    return messages_new
end



function comp_error(G::Int64, messages::Array{Float64, 1}, 
                    messages_new::Array{Float64, 1})
    max_delta = -1
    sdelta=0
        for comb_i in 1:2^G
            delta = abs(messages_new[comb_i] - messages[comb_i]) /
                            messages[comb_i]
            sdelta+=abs(delta)
            messages[comb_i] = messages_new[comb_i]
            if delta > max_delta
                max_delta = delta
            end
        end
        return max_delta, messages
end


function update_messages(temp::Float64, c::Int64, G::Int64, alpha::Float64, 
                         spins_combs::Array{Array{Int8, 1}, 1}, messages::Array{Float64, 1}, 
                         messages_new::Array{Float64, 1})
    
    messages_new = sum_prod_rule(temp, c, G, alpha, spins_combs, messages, messages_new)
    error, messages = comp_error(G, messages, messages_new)
    return error
end


# Computes the energy <s_i . s_j>
function ener(G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
              messages::Array{Float64, 1}, alpha::Float64, temp::Float64, c::Int64)
    cumul_prob = 0.0
    cumul_edge = 0.0
    for comb_i in 1:2^G
        for comb_j in 1:2^G
            d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])
            hamilt = ((2 * alpha - 1) * d_prod + abs(d_prod)) / 2 / G
            exp_hamilt = exp(hamilt / temp)
            cumul_edge += messages[comb_i] * messages[comb_j] * hamilt * exp_hamilt
            cumul_prob += messages[comb_i] * messages[comb_j] * exp_hamilt
        end
    end
    return cumul_edge / cumul_prob, cumul_prob
end


function mag_towards_one(c::Int64,G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
    messages::Array{Float64, 1}, alpha::Float64, temp::Float64)
    cumul_num = 0.0
    cumul_den = 0.0
    for comb_i in 1:2^G
        cumulant = 0.0
        for comb_j in 1:2^G
            d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])   # This is s_i . s_j
            prod = messages[comb_j] * exp(((2 * alpha - 1) * d_prod + abs(d_prod)) / temp / 2 / G) 
            cumulant += prod
        end
        probi = cumulant^c

        norm_d_prod = dot(spins_combs[comb_i], spins_combs[1]) / G
        cumul_num += norm_d_prod * probi
        cumul_den += probi
    end
    return cumul_num / cumul_den
end


function mag_norm(G::Int64, spins_combs::Array{Array{Int8, 1}, 1}, 
    messages::Array{Float64, 1}, c::Int64, alpha::Float64, temp::Float64)
    cumul_num = zeros(G)
    cumul_den = 0.0
    for comb_i in 1:2^G
        cumulant = 0.0
        for comb_j in 1:2^G
            d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])   # This is s_i . s_j
            prod = messages[comb_j] * exp(((2 * alpha - 1) * d_prod + abs(d_prod)) / temp / 2 / G) 
            cumulant += prod
        end
        probi = cumulant^c

        cumul_num .+= spins_combs[comb_i] * probi
        cumul_den += probi
    end
    av_vector = cumul_num / cumul_den

    norm = sqrt(sum(av_vector .^ 2) / G)
    return norm, av_vector, cumul_den
end



# This funtion converts an integer into an array of elements 1 and -1.
function int_to_spins(comb::Int64, G::Int64)
    str = bitstring(comb)   # The integer is encoded into a string of 0 and 1 (bits)
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


function belief_propagation(temp::Float64, c::Int64, max_iter::Int64,
    delta::Float64, G::Int64, alpha::Float64, spins_combs::Array{Array{Int8, 1}, 1}, 
    messages::Array{Float64, 1}, messages_new::Array{Float64, 1})
    
    iter = 0
    max_delta = delta + 1
    while (iter < max_iter && max_delta > delta)
        max_delta = update_messages(temp, c, G, alpha, spins_combs, messages, messages_new)
        iter += 1
    end
    return max_delta < delta
end



# This function runs several BP convergences for different parameters. In this case, the parameters 
# are T and α
function run_BP(c, seed, string1, string_init, max_iter, delta, G, 
                                temp, d_alpha, cond_init)
    alpha_min = 0.5
    alpha_max = 1.0
    alphas_list = range(alpha_min, alpha_max, step = d_alpha) 
    spins_combs = get_spins_combs(G)
    fileeq = string("BP_", string1,"_allalpha_eq_c_", string(c), "_T_", string(temp), 
                        "_delta_", string(exponent), "_maxiter_", 
                        string(max_iter), string_init, "_G_", string(G), ".txt")
    filemag_comp = string("BP_", string1,"_allalpha_mag_comp_c_", string(c), "_T_", string(temp), 
                          "_delta_", string(exponent), "_maxiter_", 
                          string(max_iter), string_init, "_G_", string(G), ".txt")
    w = open(fileeq, "w")
    f_comp = open(filemag_comp, "w") 
    write(w, "#  T   alpha   e   m(norm)    m(1)   free_ener   convergence\n")
    write(f_comp, "#  T   alpha   mag_components\n")
    for alpha in alphas_list
        if cond_init == "rand"
            messages, messages_new = set_messages_rand(seed, G)
        elseif cond_init == "ord"
            messages, messages_new = set_messages_ord(seed, G)
        end
        conv = belief_propagation(temp, c, max_iter, delta, G, alpha, spins_combs, 
                                  messages, messages_new)
        if !conv
            println("BP did not converge for T=", temp, " and alpha=",alpha)
        end
        eq_e, zij = ener(G, spins_combs, messages, alpha, temp, c)
        eq_m_norm, mag_components, zi = mag_norm(G,spins_combs,messages,c, alpha, temp)
        eq_m_oriented = mag_towards_one(c, G, spins_combs, messages, alpha, temp)
        f = -temp * (log(zi) - c * log(zij) / 2) 
        write(w, string(temp) * "\t" * string(alpha) * "\t" * string(eq_e) * "\t" *string(eq_m_norm) 
                 * "\t" *string(eq_m_oriented) * "\t" * string(f) * "\t" * string(conv) * "\n")
        write(f_comp, string(temp) * "\t" * string(alpha))
        for comp in mag_components
            write(f_comp, "\t" * string(comp))
        end
        write(f_comp, "\n")
    end
    close(w)
    close(f_comp)
end

function map_all_states(c, seed, string1, string_init, max_iter, delta, G, d_alpha, 
                               temp_list, cond_init)
    for i in eachindex(temp_list)
        temp = temp_list[i]
        run_BP(c, seed, string1, string_init, max_iter, delta, G, temp, d_alpha, cond_init)
    end
end


c = parse(Int64, ARGS[1])      # Graph connectivity

exponent = parse(Int64, ARGS[2])
delta = 10.0 ^ exponent   # When the error is below this threshold, I declare that BP converges
max_iter = 100000            # BP will stop after "max_iter" iterations
string1 = "RGfc"          # An identification of the graph you are using. This string is insterted in 
                          # the output strings
G = parse(Int64, ARGS[3])


d_alpha = 0.005

temp_list = [0.1, 0.2, 0.3, 0.4, 0.5]

cond_init = ARGS[4]
seed = parse(Int64, ARGS[5])                  # The seed for the initial messages in case 
# cond_init is 'rand', or the specific initial configuration in case cond_init=ord

println("cond_init=" * cond_init)

if cond_init == "rand"
    string_init = string("_rand_init_seed_", string(seed))
elseif cond_init == "ord"
    string_init = string("_ord_init_conf_", string(seed))
else
    println("The initial condition variable 'cond_init' must be 'ord' or 'rand'")
    exit(1)
end

map_all_states(c, seed, string1, string_init, max_iter, delta, G, d_alpha,temp_list, cond_init)
