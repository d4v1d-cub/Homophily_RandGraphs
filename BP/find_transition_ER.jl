include("BP_popdyn_ER.jl")


cm = parse(Float64, ARGS[1])      # Graph connectivity

delta = parse(Float64, ARGS[2])  # When the error is below this threshold, I declare that BP converges

if tryparse(Float64, ARGS[3]) !== nothing
    cond_init = parse(Float64, ARGS[3])
else
    cond_init = ARGS[3]
end

G = parse(Int64, ARGS[4])      # Number of opinions
S = parse(Int64, ARGS[5])      # Population size
smooth_size = parse(Int64, ARGS[6])  # Number of iterations to smooth the population dynamics
max_iter = parse(Int64, ARGS[7])            # BP will stop after "max_iter" iterations
# alpha is a parameter of the model, which is the probability of a positive interaction between two agents
# it will take values between alpha_min and alpha_max, with step alpha_step. The program stops the first time
# it finds a nonzero magnetization
alpha_min = parse(Float64, ARGS[8])
alpha_max = parse(Float64, ARGS[9])
alpha_step = parse(Float64, ARGS[10]) 
temp = parse(Float64, ARGS[11])  # Temperature of the system, which is the inverse of the noise in the model
measure_size = parse(Int64, ARGS[12])   # Number of measures to average the final results
seed = parse(Int64, ARGS[13]) # Seed for the random number generator
delta_mag = parse(Float64, ARGS[14]) # Threshold for the magnetization to be considered nonzero
path_output = ARGS[15] # Path to the output file
fileoutput = ARGS[16] # name of the output file


fileeq = string(path_output, "/", fileoutput)               
w = open(fileeq, "w")
write(w, "#  T   alpha   e   m(norm)   f    convergence\n")

for alpha in alpha_min:alpha_step:alpha_max
    println("\nRunning BP for alpha = ", alpha)
    flush(stdout)
    
    # Run the BP algorithm to find the transition point
    Random.seed!(seed)
    mag_av, ener_av, f_av = run_BP_find_transition(cm, max_iter, delta, G, temp, alpha, S, cond_init, smooth_size, measure_size)
    write(w, string(temp) * "\t" * string(alpha) * "\t" * string(ener_av) * "\t" * string(mag_av) *"\t" * string(f_av) * "\n")
    println("Final results for T=", temp, " and alpha=",alpha)
    println(string(ener_av) * "\t" * string(mag_av) *"\t" * string(f_av))
    flush(stdout)
    
    if mag_av > delta_mag
        println("Magnetization is nonzero for alpha = ", alpha)
        println("Stopping...")
        flush(stdout)
        break
    end
    println("finished alpha = ", alpha)
end

close(w)
    
