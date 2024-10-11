using LightGraphs, LinearAlgebra 
import Random
include("BP_RRG_homophily_single_eq.jl")


function map_all_states(c, param, string1, string_init, max_iter, delta, G, d_alpha, 
                               temp_list, cond_init)
    fileeq = string("BP_", string1,"_allalpha_eq_c_", string(c), 
                    "_delta_", string(exponent), "_maxiter_", 
                    string(max_iter), string_init, "_G_", string(G), ".txt")
    filemag_comp = string("BP_", string1,"_allalpha_mag_comp_c_", string(c), 
                          "_delta_", string(exponent), "_maxiter_", 
                          string(max_iter), string_init, "_G_", string(G), ".txt")
    w = open(fileeq, "w")
    f_comp = open(filemag_comp, "w")
    write(w, "#  T   alpha   e   m(norm)    m(1)   free_ener   convergence\n")
    write(f_comp, "#  T   alpha   mag_components\n")
    for temp in temp_list
        alpha_min = 0.5
        alpha_max = 1.0
        alphas_list = range(alpha_min, alpha_max, step = d_alpha) 
        spins_combs = get_spins_combs(G)
        for alpha in alphas_list
            if cond_init == "rand"
                messages, messages_new = set_messages_rand(param, G)
            elseif cond_init == "ord"
                messages, messages_new = set_messages_ord(param, G)
            elseif cond_init == "hom"
                messages, messages_new = set_messages_hom(G, param)
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
    end
    close(w)
    close(f_comp)
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

println("cond_init=" * cond_init)

if cond_init == "rand"
    param = parse(Int64, ARGS[5])                  # The seed for the initial messages in case 
                                                   # cond_init is 'rand'
    string_init = string("_rand_init_seed_", string(param))
elseif cond_init == "ord"
    param = parse(Int64, ARGS[5])                  # or the specific initial configuration in case cond_init=ord
    string_init = string("_ord_init_conf_", string(param))
elseif cond_init == "hom"
    param = parse(Float64, ARGS[5])                  # or the specific initial configuration in case cond_init=ord
    string_init = string("_hom_init_p_", string(param))
else
    println("The initial condition variable 'cond_init' must be 'ord', 'rand', or 'hom'")
    exit(1)
end

map_all_states(c, param, string1, string_init, max_iter, delta, G, d_alpha,temp_list, cond_init)
