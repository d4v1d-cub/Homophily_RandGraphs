using LightGraphs, LinearAlgebra 
import Random
include("BP_RRG_homophily_single_eq.jl")



# This function runs several BP convergences for different parameters. In this case, the parameters 
# are T and α
function run_BP_find_transition(c, seed, string1, string_init, max_iter, delta, G, 
                                temp, d_alpha, m_thr, cond_init)
    alpha_min = 0.5
    alpha_max = 1.0
    alphas_list = range(alpha_min, alpha_max, step = d_alpha) 
    spins_combs = get_spins_combs(G)
    fileeq = string("BP_", string1,"_eq_c_", string(c), "_T_", string(temp), 
                        "_delta_", string(exponent), "_maxiter_", 
                        string(max_iter), string_init, "_G_", string(G), ".txt")
    filemag_comp = string("BP_", string1,"_mag_comp_c_", string(c), "_T_", string(temp), 
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
            messages, messages_new = set_messages_hom(seed, G)
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
        
        if abs(eq_m_norm) > m_thr && conv
            return alpha
            break
        end
    end
    println("No transition found in the interval [", alpha_min, ", ", alpha_max, "]")
    close(w)
    close(f_comp)
    return 1
end

function numerical_transitions(c, seed, string1, string_init, max_iter, delta, G, d_alpha, m_thr, 
                               fileout,temp_list, cond_init)
    w = open(fileout, "w")
    for i in eachindex(temp_list)
        temp = temp_list[i]
        alpha_num = run_BP_find_transition(c, seed, string1, string_init, max_iter, 
                                           delta, G, temp, d_alpha, m_thr, cond_init)
        println("\n" * "Numerical transition:  T=" * string(temp) * "\talpha=" * string(alpha_num) * "\n")
        write(w, string(temp) * "\t" * string(alpha_num) * "\n")
    end
    close(w)
end


c = parse(Int64, ARGS[1])      # Graph connectivity

exponent = parse(Int64, ARGS[2])
delta = 10.0 ^ exponent   # When the error is below this threshold, I declare that BP converges
max_iter = 100000            # BP will stop after "max_iter" iterations
string1 = "RGfc"          # An identification of the graph you are using. This string is insterted in 
                          # the output strings
G = parse(Int64, ARGS[3])


d_alpha = 0.001
m_thr = 0.01

temp_range=range(0.01,stop=0.54,step=0.01)
temp_list = collect(temp_range)

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

fileout = string("Phase_Diagram_Homophily_numerical_c_", string(c), 
            "_delta_", string(exponent), "_maxiter_", 
            string(max_iter), string_init, "_G_", string(G), ".txt")

numerical_transitions(c, seed, string1, string_init, max_iter, delta, G, d_alpha, m_thr, fileout,temp_list,
                      cond_init)
