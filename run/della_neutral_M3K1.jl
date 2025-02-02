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
const M = 3 # number of sets
const K = 1 # number of issues to care about
const P = 2 # number of parties

const β = 0.000    # selection strength
# const u = 0.001    # strategy mutation rate
# const v = 0.001    # set mutation rate
# const p = 0.0     # party bias (0 = low, 1 = high)
# const ϵ = 1.       # damping factor (0 = no bias in mutation, 1 = max bias (= p) in mutation)
# note: α = 0 automatically

const pvals = repeat([0.0, 0.25, 0.5, 0.75, 1.0], 1) # uniform populations for now
const us    = [0.001]
const vs    = [0.001, 10^((log10(0.001)+log10(0.005))/2), 0.005, 10^((log10(0.005)+log10(0.025))/2), 0.025, 10^((log10(0.025)+log10(0.125))/2), 0.125, 10^((log10(0.125)+log10(0.625))/2), 0.625] 
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
total_summary = SharedArray{Float64, 2}(length(pvals)*length(us)*length(vs),26) # SharedArray{Real}[] 

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
    n_tracker = NeutralTracker(pop, generations, verbose_track)

    println("\nstarting $(filename)")

    ##########################
    ## main simulation loop ##
    ##########################
    for gen::Int64 in 1:generations
        # if generations > 10 && (g-1) % (generations ÷ 10) == 0 println("\tinitiating generation $g") end
        # update population
        evolve_neutral!(n_tracker.pop, 1) # evolve without selection
        track_neutral!(n_tracker, gen)

        if gen == generations 
            println("replicate complete, ran $(generations) generations:\n") 
            println("\tnum_imitating:\t\t\t$(n_tracker.pop.num_imitating)")
            println("\tnum_not_imitating:\t\t$(n_tracker.pop.num_not_imitating)")
            println("\tnum_strategies_changed:  \t$(n_tracker.pop.num_strategies_changed)")
            println("\t\tnum_random_strategies: \t$(n_tracker.pop.num_random_strategies)")
            println("\tnum_sets_changed: \t\t$(n_tracker.pop.num_sets_changed)")
            println("\t\tnum_random_opinions: \t$(n_tracker.pop.num_random_opinions)")
            println("\t\tnum_biased_opinions: \t$(n_tracker.pop.num_biased_opinions)")
            println("\tnum_both_changed: \t\t$(n_tracker.pop.num_both_changed)\n")
        end
    end

    ############################
    ## prep data for plotting ##
    ############################
    # make a summary table
    # ens_summary, colnames = make_summary_table( tracker, gen_cutoff, generations )
    # total_summary[i,:]    = ens_summary

    ens_summary = hcat( generations, sets.N, sets.M, sets.K, sets.P,                            # population parameters
                        game.b, game.c, game.β, game.u, game.v, game.ps[1], game.ps[2], game.ϵ, # game parameters
                        n_tracker.ens_strat_running[generations,:]', # strategy distribution
                        n_tracker.ens_y_running[generations],
                        n_tracker.ens_z_running[generations],
                        n_tracker.ens_g_running[generations],
                        n_tracker.ens_h_running[generations],
                        n_tracker.ens_zp_running[generations],
                        n_tracker.ens_gp_running[generations],
                        n_tracker.ens_hp_running[generations],
                        n_tracker.ens_sia_sid_running[generations],
                        n_tracker.ens_sia_sjd_running[generations]
                        )
    col_list = hcat("generations", "N", "M", "K", "P",        # population parameters
                    "b", "c", "β", "u", "v", "p1", "p2", "ϵ", # game parameters
                    "DD_mean_final", "DC_mean_final", "CD_mean_final", "CC_mean_final",
                    "y", "z", "g", "h", "zp", "gp", "hp", "sia_sid", "sia_sjd"
                    )
    colnames = join(col_list, ',')
    total_summary[i,:] = ens_summary

    # overwrite file after every parameter set (in case code stops due to memory issues)
    # writedlm("data/$(testtype)_run1.csv"; [colnames; total_summary], ',', header=true)
    open("data/gens_$(generations)/$(testtype)_neutral_M$(M)K$(K)_$(start_time).csv", "w") do f
        write(f, "$(colnames)\n")
        writedlm(f, total_summary, ',')
    end
    println("data processed and exported")

end
