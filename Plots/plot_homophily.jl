using Plots
using LaTeXStrings

function ImportData(filename)
    r = open(filename, "r")
    vals_alpha = Array{Float64, 1}()
    vals_e = Array{Float64, 1}()
    vals_rho = Array{Float64, 1}()
    while !eof(r)
        line = split(readline(r))
        if line[end] == "true"
            push!(vals_alpha, parse(Float64, line[2]))
            push!(vals_e, parse(Float64, line[3]))
            push!(vals_rho, parse(Float64, line[4]))
        end
    end
    close(r)
    return vals_alpha, vals_e, vals_rho
end


function ImportAll(direc, head_str, common_str_1, common_str_2, temps)
    all_alpha = Array{Array{Float64, 1}, 1}()
    all_e = Array{Array{Float64, 1}, 1}()
    all_rho = Array{Array{Float64, 1}, 1}()
    for temp in temps
        filename = direc * head_str * "_" * common_str_1 * "_T_" * string(temp) * 
                   common_str_2 * ".txt"
        vals_alpha, vals_e, vals_rho = ImportData(filename)
        push!(all_alpha, vals_alpha)
        push!(all_e, vals_e)
        push!(all_rho, vals_rho)
    end
    return all_alpha, all_e, all_rho
end


function plot_several_temps(x, y, filefig, temps, ylabel)
    temp = temps[1]
    cl = color_list(:tab10)
    scatter(x[1], y[1], label="T = " * string(temp), lw=3, mc=cl[1])
    plot!(x[1], y[1], label="", lw=1, lc=cl[1])

    for i in 2:length(x)
        scatter!(x[i], y[i], label="T = " * string(temps[i]), lw=3, mc=cl[i])
        plot!(x[i], y[i], label="", lw=1, lc=cl[i])
    end

    title!("G=" * string(G))
    xlabel!(L"\alpha")
    ylabel!(ylabel)
    savefig(filefig)

end


function plot_all(temps, direc, head_str, common_str_1, common_str_2)
    
    all_alpha, all_e, all_rho = ImportAll(direc, head_str, common_str_1, common_str_2, temps)

    filefig = head_str * "_ener_" * common_str_1 * common_str_2 * ".png"
    plot_several_temps(all_alpha, all_e, filefig, temps, "e")
    
    filefig = head_str * "_rho_pp_" * common_str_1 * common_str_2 * ".png"
    plot_several_temps(all_alpha, all_rho, filefig, temps, L"\rho_{++}")
end

n = 1000
c = 3
temps = [0.5, 0.2]
exp_delta = -15
sampling = 10
max_iters = 1000
seed = 1
G = 2

direc = "Data_single_instance/"
head_str = "BP_RGfc_eq"

common_str_1 = "N_" * string(n) * "_c_" * string(c)
common_str_2 = "_delta_" * string(exp_delta) * "_sampling_" * string(sampling) * 
"_maxiter_" * string(max_iters) * "_seed_" * string(seed) * 
"_G_" * string(G)

plot_all(temps, direc, head_str, common_str_1, common_str_2)