using Interpolations


function read_pairs(fileinput::String)
    temps_data = Float64[]
    alphas_data = Float64[]
    max_ndigits = 0
    open(fileinput, "r") do f
        for line in eachline(f)
            if !startswith(line, "#") && !isempty(line)
                parts = split(line)
                push!(temps_data, parse(Float64, parts[1]))
                push!(alphas_data, parse(Float64, parts[2]))
                # Check the number of digits in the alpha value
                splitted_number = split(parts[2], '.')
                if length(splitted_number) < 2
                    continue  # No decimal point, skip this line
                else
                    ndigits = length(split(parts[2], '.')[2])  # Count digits after the decimal point
                    if ndigits > max_ndigits
                        max_ndigits = ndigits
                    end
                end
            end
        end
    end
    return temps_data, alphas_data, max_ndigits
end


function interpolate_pairs(temps_data::Array{Float64, 1}, alphas_data::Array{Float64, 1}, temp0::Float64, tempf::Float64, dtemp::Float64)
    temps = temp0:dtemp:tempf
    alphas = Float64[]
    interp_linear_extrap = linear_interpolation(temps_data, alphas_data, extrapolation_bc=Line());
    for temp in temps
        push!(alphas, interp_linear_extrap(temp))
    end
    return temps, alphas
end


function print_results(fileoutput::String, temps, alphas, ndigits::Int64)
    open(fileoutput, "w") do f
        for i in eachindex(temps)
            write(f, string(temps[i]), "\t", string(round(alphas[i]; digits=ndigits)), "\n")
        end
    end
end


fileinput = ARGS[1]  # Input file with parameters proposed parameters
fileoutput = ARGS[2] # Output file for the results
temp0 = parse(Float64, ARGS[3])  # Initial temperature
tempf = parse(Float64, ARGS[4])  # Final temperature
dtemp = parse(Float64, ARGS[5])  # Temperature step


temps_data, alphas_data, max_ndigits = read_pairs(fileinput)
temps, alphas = interpolate_pairs(temps_data, alphas_data, temp0, tempf, dtemp)
print_results(fileoutput, temps, alphas, max_ndigits)
