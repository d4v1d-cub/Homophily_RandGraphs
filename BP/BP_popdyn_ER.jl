using LightGraphs, LinearAlgebra
using Random
using Distributions
using Statistics



function int_to_spins(comb::Int64, G::Int64)
    str = bitstring(comb)   
    spins = Array{Int8, 1}()
    for i in 0:G-1
        push!(spins, 1 - 2 * parse(Int64, str[end-i]))   
    end
    return spins
end

# Compute all possible 1 and -1 combinations given a certain dimension

function get_spins_combs(G::Int64)
    spins_combs = Array{Array{Int8, 1}, 1}()
    for comb in 1:2^G
        push!(spins_combs, int_to_spins(comb - 1, G))
    end
    return spins_combs
end

# Normalize a vector

function normalize_norm_1(vec::Array{Float64, 1})
    norm = sum(vec)
    for i in 1:length(vec)
        vec[i] /= norm
    end
    return vec
end


# Generate a population of messages, each of them with random components,
# each message is also normalized

function set_population_rand(G::Int64,S::Int64)
    pop=Array{Array{Float64, 1}, 1}()
    for _=1:S
        push!(pop, normalize_norm_1([rand() for j=1:2^G]))        
    end
    return pop
end

function set_population_ord(G::Int64,S::Int64)
    pop=Array{Array{Float64, 1}, 1}()
    for i=1:S
        push!(pop, zeros(2^G))
        pop[i][1]=1
        pop[i]=normalize_norm_1(pop[i])        
    end
    return pop
end

function set_population_arbitrary(G::Int64,S::Int64,p::Float64)
    pop=Array{Array{Float64, 1}, 1}()
    for i=1:S
        messages = Array{Float64, 1}()
        mess = Array{Float64, 1}()
        combs=get_spins_combs(G)
        for i in 1:2^G
            kexp=sum(x for x in combs[i] if x < 0)
            bin_p=p^(G-abs(kexp))*(1-p)^abs(kexp)
            push!(mess, bin_p)
        end
   	    messages=deepcopy(mess)
        push!(pop, messages)
    end
    return pop
end

# Computes a new message given a number c number of messages

function sum_prod_rule(ncombs::Int64, random_mess::Array{Array{Float64, 1}, 1}, exp_factor::Array{Float64, 2})
    message=zeros(ncombs)    
    for comb_i in 1:ncombs
        prod = 1.0
        for r_mess in random_mess
            cumulant = 0.0
            for comb_j in 1:ncombs
                cumulant += r_mess[comb_j] * exp_factor[comb_i, comb_j]  # ηᵏ(σₖ) * ψₖ(σ, σₖ) 
            end
            prod *= cumulant
        end
        message[comb_i] = prod
    end
    return normalize_norm_1(message)
end


# Computes the energy <s_i . s_j>
function ener(ncombs::Int64, mij::Array{Float64, 1}, mji::Array{Float64, 1}, hamilt::Array{Float64, 2}, exp_factor::Array{Float64, 2})
    cumul_prob = 0.0
    cumul_edge = 0.0
    for comb_i in 1:ncombs
        for comb_j in 1:ncombs
            cumul_edge += mij[comb_i] * mji[comb_j] * hamilt[comb_i, comb_j] * exp_factor[comb_i, comb_j]
            cumul_prob += mij[comb_i] * mji[comb_j] * exp_factor[comb_i, comb_j]
        end
    end
    return cumul_edge / cumul_prob, cumul_prob
end



function ener_eval_final(pop::Array{Array{Float64, 1}, 1}, ncombs::Int64, hamilt::Array{Float64, 2}, exp_factor::Array{Float64, 2})
    cumul_log_zij = 0
    cumul = 0
    for _=1:length(pop)
        rand_mess = Distributions.sample(pop, 2)
        eq_e, zij = ener(ncombs, rand_mess[1], rand_mess[2], hamilt, exp_factor)
        cumul += eq_e
        cumul_log_zij += log(zij)
    end
    energy = cumul / length(pop)
    log_z_ij = cumul_log_zij / length(pop)
    
    return energy, log_z_ij 
end



function mag_norm(ncombs::Int64, spins_combs::Array{Array{Int8, 1}, 1}, rand_mess::Array{Array{Float64, 1}, 1}, exp_factor::Array{Float64, 2})

    cumul_num = zeros(G)
    cumul_den = 0.0
    for comb_i in 1:ncombs  
        prod = 1.0
        for r_mess in rand_mess  
            cumulant = 0.0
            for comb_j in 1:ncombs
                cumulant += r_mess[comb_j] * exp_factor[comb_i, comb_j]  # ηᵏ(σₖ) * ψₖ(σ, σₖ)
            end
            prod *= cumulant 
        end

        cumul_den += prod
        cumul_num .+= spins_combs[comb_i] * prod
    end
    av_vector = cumul_num / cumul_den

    norm = sqrt(sum(av_vector .^ 2) / G)
    return norm, cumul_den
end



function mag_pop_dynam(cm::Float64, pop::Array{Array{Float64, 1}, 1}, ncombs::Int64, spins_combs::Array{Array{Int8, 1}, 1}, exp_factor::Array{Float64, 2})
    cumul = 0
    for _ in 1:length(pop)
        c = rand(Poisson(cm))
        rand_mess = Distributions.sample(pop, c)
        eq_m_norm, _ = mag_norm(ncombs, spins_combs, rand_mess, exp_factor)
        cumul += eq_m_norm
    end
    return cumul/length(pop)
end


function mag_eval_final(cm::Float64, pop::Array{Array{Float64, 1}, 1}, ncombs::Int64, spins_combs::Array{Array{Int8, 1}, 1}, exp_factor::Array{Float64, 2})
    cumul = 0
    cumul_log_zi = 0
    for _=1:length(pop)
        c = rand(Poisson(cm))
        rand_mess = Distributions.sample(pop, c)
        eq_m_norm, zi = mag_norm(ncombs, spins_combs, rand_mess, exp_factor)
        cumul += eq_m_norm
        cumul_log_zi += log(zi)
    end
    mag = cumul / length(pop)
    log_z_i = cumul_log_zi / length(pop)
    
    return mag, log_z_i
end


function popdyn(S::Int64, cm::Float64, alpha::Float64, delta::Float64, ncombs::Int64, spins_combs::Array{Array{Int8, 1}, 1}, pop::Array{Array{Float64, 1}, 1}, 
                max_iter::Int64, smooth_size::Int64, exp_factor::Array{Float64, 2}; conv_param::Int64=5, iters_print::Int64=2, bool_print::Bool=true)

    time_last_print = time()
    av_time_iter = 0.0
    iter = 0
    mag_av_prev = 0.0   
    convergence_counter = 0
    
    while (iter < max_iter && convergence_counter < conv_param)
        start_time = time()
        mag_av = 0.0
        for _ in 1:smooth_size
            for _ in 1:S
                c = rand(Poisson(cm))
                rand_mess = Distributions.sample(pop, c)
                r = rand(1:length(pop))
                pop[r] = sum_prod_rule(ncombs, rand_mess, exp_factor)
            end
            single_mag = mag_pop_dynam(cm, pop, ncombs, spins_combs, exp_factor)
            mag_av += single_mag
        end
        mag_av /= smooth_size
        dmag = abs(mag_av - mag_av_prev)
        mag_av_prev = mag_av
                
        if (dmag < delta)
            convergence_counter += 1
        else
            convergence_counter = 0
        end
                
        
        av_time_iter += time() - start_time
        iter += 1
        
        if iter % iters_print == 0 && bool_print
            time_between_prints = time() - time_last_print    
            time_last_print = time()
            av_time_iter /= iters_print
            println("Average time last $iters_print iterations: $av_time_iter  seconds")
            println("Total Time for last $iters_print iterations: $time_between_prints  seconds")
            println("mag = $mag_av   dmag = $dmag")
            println("alpha = $alpha")
            flush(stdout)
            av_time_iter = 0.0 
        end
    end
    conv = (convergence_counter >= conv_param)
    if bool_print
        println("Converge for $alpha:  $conv")
        flush(stdout)
    end
    return conv, pop
end


function popdyn_measure(S::Int64, cm::Float64, ncombs::Int64, spins_combs::Array{Array{Int8, 1}, 1}, pop::Array{Array{Float64, 1}, 1}, 
                        measure_size::Int64, hamilt::Array{Float64, 2}, exp_factor::Array{Float64, 2})
     
    mag_av = 0.0
    ener_av = 0.0
    f_av = 0.0
    for _ in 1:measure_size
        for _ in 1:S
            c = rand(Poisson(cm))
            rand_mess = Distributions.sample(pop, c)
            r = rand(1:length(pop))
            pop[r] = sum_prod_rule(ncombs, rand_mess, exp_factor)
        end
        single_mag, log_zi = mag_eval_final(cm, pop, ncombs, spins_combs, exp_factor)
        single_ener, log_zij = ener_eval_final(pop, ncombs, hamilt, exp_factor)  
        mag_av += single_mag
        ener_av += single_ener
        f_av += -log_zi + cm * log_zij / 2
    end
    mag_av /= measure_size
    ener_av /= measure_size
    f_av /= measure_size
    return mag_av, ener_av, f_av
end


function build_auxiliary_matrices(ncombs::Int64, G::Int64, alpha::Float64, temp::Float64, spins_combs::Array{Array{Int8, 1}, 1})
    hamilt = zeros(ncombs, ncombs)
    exp_factor = zeros(ncombs, ncombs)
    
    for comb_i in 1:ncombs
        for comb_j in 1:ncombs
            d_prod = dot(spins_combs[comb_i], spins_combs[comb_j])
            hamilt[comb_i, comb_j] = ((2 * alpha - 1) * d_prod + abs(d_prod)) / 2 / G
            exp_factor[comb_i, comb_j] = exp(hamilt[comb_i, comb_j] / temp)
        end
    end
    return hamilt, exp_factor
end


function run_BP_find_transition(cm::Float64, max_iter::Int64, delta::Float64, G::Int64, temp::Float64, alpha::Float64, S::Int64, 
                                cond_init::Union{String, Float64}, smooth_size::Int64, measure_size::Int64; bool_print::Bool=true)  
    ncombs = 2^G
    spins_combs = get_spins_combs(G)
       
    if cond_init == "rand"
        population_mess = set_population_rand(G,S)
    elseif cond_init == "ord"
        population_mess = set_population_ord(G,S)
	else
	    population_mess = set_population_arbitrary(G, S, cond_init)
    end
    
    
    hamilt, exp_factor = build_auxiliary_matrices(ncombs, G, alpha, temp, spins_combs)
    conv, population_mess = popdyn(S, cm, alpha, delta, ncombs, spins_combs, population_mess, max_iter, smooth_size, exp_factor, bool_print=bool_print)
        
    if !conv
        println("BP did not converge for T=", temp, " and alpha=",alpha)
        flush(stdout)
    end
        
    mag_av, ener_av, f_av = popdyn_measure(S, cm, ncombs, spins_combs, population_mess, measure_size, hamilt, exp_factor)
    
    return mag_av, ener_av, f_av
end

    
