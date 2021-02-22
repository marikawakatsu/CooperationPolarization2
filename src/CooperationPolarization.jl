# __precompile__()
module CooperationPolarization

# structs in setup.jl
export Sets, Game, Simulation, Population
export Tracker

# functions in setup.jl
export compute_total_memberships
export party_assignment
export set_members_and_pairs
export random_sets
export vectobase3

# functions in simulate.jl
export evolve!
export update_payoffs_using_matrix!
export update_strategies_and_opinions_db!
export update_opinions!
export update_strategies!
export no_strategy_or_opinion_update!
export compute_payoffs!

# functions in plot.jl
export plot_strat_dist
export plot_opn_dist
export plot_set_dist
export plot_total_opn
export base3tovec
export plotter

# functions in track.jl
export track!
export track_payoff_means!
export track_strat_dist!
export track_pop_dist!
export track_running_strat_dist!
export track_opn_dist!
export track_total_means_and_vars!
export track_total_cooperation!
export track_opn_simpson!
export track_sets_simpson!
export track_cityblock!
export track_hamming!
export make_summary_table
export make_col_list
export make_ens_summary

# functions in simulate_neutral.jl
export evolve_neutral!
export update_actions_neutral!
export compute_actions_neutral!
export update_strategies_and_opinions_neutral!

# struct and functions in track_neutral.jl
export NeutralTracker
export track_neutral!

# load 
include("setup.jl")
include("simulate.jl")
include("track.jl")
include("plot.jl")
include("simulate_neutral.jl")
include("track_neutral.jl")

end
