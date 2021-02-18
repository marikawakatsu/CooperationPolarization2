# sim_multithread_della.jl: main simulation code with multithreading
# on della
using Distributed
using Dates
using Pkg; Pkg.activate(".")
using CooperationPolarization
addprocs(15; exeflags="--project")

@everywhere begin
    ###################
    ## load packages ##
    ###################
    using Distributed
    using Dates
    using Pkg; Pkg.activate(".")
    using CooperationPolarization
    using BenchmarkTools
    using PyPlot: plt, savefig
    using Statistics: mean, var
    using SharedArrays
    using DelimitedFiles
    println("packages loaded")

end

###########################
## initialize parameters ##
###########################
# model parameters (fixed)
const b = 1.0 # benefit to cooperating
const c = 0.2 # cost to cooperating

const N = 40 # population size
const P = 2 # number of parties

const β = 0.001    # selection strength
# const u = 0.001    # strategy mutation rate
# const v = 0.001    # set mutation rate
# const p = 0.       # party bias (0 = low, 1 = high)
# const ϵ = 1.       # damping factor (0 = no bias in mutation, 1 = max bias (= p) in mutation)
# note: α = 0 automatically

# model parameters (varied)
const M = 5 # number of sets
const K = 4 # number of issues to care about

const pvals = repeat([0.0, 0.25, 0.5, 0.75, 1.0], 20) # uniform populations for now
const us    = [0.001]
const vs    = [0.005]
const ϵ     = 1.
const parameters = collect(Iterators.product(pvals,us,vs))

# simulation parameters
const verbose       = false
const verbose_track = false
const generations   = 20000000
const testtype      = "della"  # specify location

# plotting and data processing parameters
const gen_skip   = 1 # max( div(generations, 10000), 1)
const gen_start  = max( 1, generations-500000 )
const gen_end    = generations
const gen_cutoff = 2000000

println("parameters defined")

const start_time = Dates.format(now(),"yymmdd-HHMMSS.sss")

################
## simulation ##
################
# initialize global tracker
total_summary = SharedArray{Float64, 2}(length(pvals)*length(us)*length(vs),543) # SharedArray{Real}[] 

@time @sync @distributed for i in 1:length(pvals)*length(us)*length(vs)

    # set parameter(s)
    p::Float64 = parameters[i][1]
    u::Float64 = parameters[i][2]
    v::Float64 = parameters[i][3]

    # set file name and type
    filename = "b_$(b)_c_$(c)_N_$(N)_M_$(M)_K_$(K)_P_$(P)_beta_$(β)_u_$(u)_v_$(v)_p_$(p)_ϵ_$(ϵ)_gens_$(generations)_$(start_time)"
    
    ########################################
    ## initialize population and trackers ##
    ########################################
    sets    = random_sets(N, M, K, P)
    game    = Game(b, c, β, u, v, p, ϵ)
    pop     = Population(sets, game, verbose)
    tracker = Tracker(pop, generations, verbose_track)

    println("\nstarting $(filename)")

    ##########################
    ## main simulation loop ##
    ##########################
    for g::Int64 in 1:generations
        # if generations > 10 && (g-1) % (generations ÷ 10) == 0 println("\tinitiating generation $g") end

        evolve!(tracker.pop, 1)
        track!(tracker, g)

        if g == generations 
            println("replicate complete, ran $(generations) generations:\n") 
            println("\tnum_imitating:\t\t\t$(tracker.pop.num_imitating)")
            println("\tnum_not_imitating:\t\t$(tracker.pop.num_not_imitating)")
            println("\tnum_strategies_changed:  \t$(tracker.pop.num_strategies_changed)")
            println("\t\tnum_random_strategies: \t$(tracker.pop.num_random_strategies)")
            println("\tnum_sets_changed: \t\t$(tracker.pop.num_sets_changed)")
            println("\t\tnum_random_opinions: \t$(tracker.pop.num_random_opinions)")
            println("\t\tnum_biased_opinions: \t$(tracker.pop.num_biased_opinions)")
            println("\tnum_both_changed: \t\t$(tracker.pop.num_both_changed)\n")
        end
    end

    ############################
    ## prep data for plotting ##
    ############################
    # make a summary table
    ens_summary, colnames = make_summary_table( tracker, gen_cutoff, generations )
    total_summary[i,:]    = ens_summary

    # overwrite file after every parameter set (in case code stops due to memory issues)
    # writedlm("data/$(testtype)_run1.csv"; [colnames; total_summary], ',', header=true)
    open("data/gens_$(generations)/$(testtype)_run_multi_M$(M)K$(K)_$(start_time).csv", "w") do f
        write(f, "$(colnames)\n")
        writedlm(f, total_summary, ',')
    end
    println("data processed and exported")

end

#################################
## plotting simulation results ##
#################################
# # plotting
# fig = plotter( tracker, gen_start, gen_skip, gen_end );

# plt.suptitle("Parameters: b=$b, c=$c, N=$N, M=$M, K=$K, β=$β, u=$u, v=$v, q=$q, gamma=$gamma", fontsize=12);
# fig.tight_layout(rect=[0, 0.03, 1, 0.97]);

# # save file
# savefig("plots/gens_$(generations)/$(testtype)_$(filename).pdf")
# plt.close(fig);
# println("plots saved")
