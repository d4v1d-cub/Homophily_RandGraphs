include("BP_popdyn_ER.jl")


cm = parse(Float64, ARGS[1])      # Graph connectivity

delta = parse(Float64, ARGS[2])  # When the error is below this threshold, I declare that BP converges

if tryparse(Float64, ARGS[3]) !== nothing
    string_init="ord_p_" * ARGS[3]
    cond_init = parse(Float64, ARGS[3])
else
    cond_init = ARGS[3]
    if cond_init == "rand"
        string_init = "_rand_init"
    elseif cond_init == "ord"
        string_init = "_init_ord"
    else
        println("The initial condition variable 'cond_init' must be 'ord', 'rand' or a floating point number")
        exit(1)
    end
end

G = parse(Int64, ARGS[4])      # Number of opinions
S = parse(Int64, ARGS[5])      # Population size
smooth_size = parse(Int64, ARGS[6])  # Number of iterations to smooth the population dynamics
max_iter = parse(Int64, ARGS[7])            # BP will stop after "max_iter" iterations
string1 = "Erdos_R"          # An identification of the graph you are using. This string is insterted in the output strings
# alpha is a parameter of the model, which is the probability of a positive interaction between two agents
alpha = parse(Float64, ARGS[8])
temp = parse(Float64, ARGS[9])  # Temperature of the system, which is the inverse of the noise in the model
measure_size = parse(Int64, ARGS[10])   # Number of measures to average the final results
seed = parse(Int64, ARGS[11]) # Seed for the random number generator
delta_mag = parse(Float64, ARGS[12]) # Threshold for the magnetization to be considered nonzero
path_output = ARGS[13] # Path to the output file
fileoutput = ARGS[14] # name of the output file
                     
fileeq = string(path_output, "/", fileoutput) 
w = open(fileeq, "w")

println("\nRunning BP for alpha = ", alpha)
flush(stdout)
    
# Run the BP algorithm to find the transition point
Random.seed!(seed)
mag_av, ener_av, f_av = run_BP_find_transition(cm, max_iter, delta, G, temp, alpha, S, cond_init, smooth_size, measure_size)
write(w, string(temp) * "\t" * string(alpha) * "\t" * string(ener_av) * "\t" * string(mag_av) *"\t" * string(f_av) * "\n")
println("Final results for T=", temp, " and alpha=",alpha)
println(string(ener_av) * "\t" * string(mag_av) *"\t" * string(f_av))
flush(stdout)
    
println("finished alpha = ", alpha)

close(w)
    
