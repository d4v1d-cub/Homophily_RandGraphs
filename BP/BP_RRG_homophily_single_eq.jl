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


function set_messages_hom(G::Int64, p::Float64)
    messages=zeros(2^G)
    for conf in 1:2^G
        bstr = bitstring(conf - 1)
        num_1 = count(x -> x == '1', bstr)
        messages[conf] = p^num_1 * (1 - p)^(G - num_1)
    end
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
