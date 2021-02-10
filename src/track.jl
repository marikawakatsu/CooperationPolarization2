# track.jl: struct and function to track population dynamics
using Statistics: sum, var, mean
using StatsBase
using Distances: pairwise, Cityblock, Hamming, cityblock, hamming

export Tracker
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

mutable struct Tracker
    # structure for storing simulation data
    pop::Population
    generations::Int64
    verbose_track::Bool                       # print messages as data are tracked
    ens_payoff_means::Array{Float64, 1}       # track average payoff in population
    ens_strat_dist::Array{Float64, 2}         # track frequencies of strategies
    ens_pop_dist::Array{Int64, 2}             # track number of individuals in each set
    ens_running_strat_dist::Array{Float64, 2} # track running average of frequencies of strategies
    ens_strat_payoff_means::Array{Float64, 2} # track average payoff for each strategy
    ens_hi::Array{Float64, 2}                 # track average opinion / set membership vector
    ens_opn_dist::Array{Float64, 2}           # track frequencies of opinions
    ens_total_means::Array{Float64, 1}        # track average of total opinions
    ens_total_vars::Array{Float64, 1}         # track variance of total opinions
    ens_total_cooperation::Array{Float64, 1}  # track effective / total cooperation
    ens_party_cooperation::Array{Float64, 1}  # track effective / total cooperation within parties
    ens_enemy_cooperation::Array{Float64, 1}  # track effective / total cooperation between parties
    ens_opn_simpson::Array{Float64,1}         # track diversity of opinions via simpson index
    ens_sets_simpson::Array{Float64,1}        # track diversity of set memberships via simpson index
    ens_cityblock_ind::Array{Float64, 1}      # track average cityblock distance between individual opinions
    ens_cityblock_party::Array{Float64, 1}    # track average cityblock distance between individuals of the same party
    ens_cityblock_enemy::Array{Float64, 1}    # track average cityblock distance between individuals of different parties
    ens_hamming_ind::Array{Float64, 1}        # track average hamming distance between individual set memberships
    ens_hamming_party::Array{Float64, 1}      # track average hamming distance between individuals of the same party
    ens_hamming_enemy::Array{Float64, 1}      # track average hamming distance between individuals of different parties

    # constructor
    function Tracker( pop::Population, generations::Int64, verbose_track::Bool )
        # initialize trackers
        ens_payoff_means       = Float64[0.0 for gens=1:generations]
        ens_strat_dist         = Float64[0.0 for gens=1:generations, strats=1:length(pop.all_strategies)]
        ens_pop_dist           = Int64[  0.0 for gens=1:generations, sets=1:pop.sets.M]
        ens_running_strat_dist = Float64[0.0 for gens=1:generations, strats=1:length(pop.all_strategies)]
        ens_strat_payoff_means = Float64[0.0 for gens=1:generations, strats=1:length(pop.all_strategies)]
        ens_hi                 = Float64[0.0 for gens=1:generations, sets=1:pop.sets.M]
        ens_opn_dist           = Float64[0.0 for gens=1:generations, opinions=1:length(pop.all_opinions)]
        ens_total_means        = Float64[0.0 for gens=1:generations]
        ens_total_vars         = Float64[0.0 for gens=1:generations]
        ens_total_cooperation  = Float64[0.0 for gens=1:generations]
        ens_party_cooperation  = Float64[0.0 for gens=1:generations]
        ens_enemy_cooperation  = Float64[0.0 for gens=1:generations]
        ens_opn_simpson        = Float64[0.0 for gens=1:generations]
        ens_sets_simpson       = Float64[0.0 for gens=1:generations]
        ens_cityblock_ind      = Float64[0.0 for gens=1:generations]
        ens_cityblock_party    = Float64[0.0 for gens=1:generations]
        ens_cityblock_enemy    = Float64[0.0 for gens=1:generations]
        ens_hamming_ind        = Float64[0.0 for gens=1:generations]
        ens_hamming_party      = Float64[0.0 for gens=1:generations]
        ens_hamming_enemy      = Float64[0.0 for gens=1:generations]

        return new(pop, generations, verbose_track,
                   ens_payoff_means, ens_strat_dist, ens_pop_dist, 
                   ens_running_strat_dist, ens_strat_payoff_means, ens_hi, 
                   ens_opn_dist, ens_total_means, ens_total_vars, 
                   ens_total_cooperation, ens_party_cooperation, ens_enemy_cooperation,
                   ens_opn_simpson, ens_sets_simpson,
                   ens_cityblock_ind, ens_cityblock_party, ens_cityblock_enemy,
                   ens_hamming_ind, ens_hamming_party, ens_hamming_enemy
                   )
    end
end

function track!( tracker::Tracker, g::Int64 )
    # main function to track the population dynamics

    if tracker.verbose_track println("generation $g") end

    track_payoff_means!(tracker, g) # track mean payoffs
    track_strat_dist!(tracker, g) # track strategy distribution
    track_pop_dist!(tracker, g) # track set distribution
    track_running_strat_dist!(tracker, g) # track running average
    track_opn_dist!(tracker, g) # track opinion distribution
    track_total_means_and_vars!(tracker, g) # track mean and variance of total opinion
    track_total_cooperation!(tracker, g, true) # track effective cooperation
    # track_opn_simpson!(tracker, g) # track simpson's index for opinions
    # track_sets_simpson!(tracker, g) # track simpson's index for set memberships
    track_cityblock!(tracker, g) # track cityblock distances
    track_hamming!(tracker, g) # track hamming distances

end

function track_payoff_means!( tracker::Tracker, g::Int64 )
    # compute and store mean payoff
    tracker.ens_payoff_means[g] = mean( tracker.pop.payoffs )     
    if tracker.verbose_track println("\t $(tracker.ens_payoff_means[g])") end
end

function track_strat_dist!( tracker::Tracker, g::Int64 )
    # compute and store frequency of individuals with each strategy
    # strategies_base10           = [2*tracker.pop.strategies[i,1] + tracker.pop.strategies[i,2] for i in 1:tracker.pop.sets.N]
    tracker.ens_strat_dist[g,:] = [count(x->x==strat, tracker.pop.strategies_base10) / tracker.pop.sets.N for strat in tracker.pop.all_strategies] 
    if tracker.verbose_track println("\t $(tracker.ens_strat_dist[g,:])") end
end

function track_pop_dist!( tracker::Tracker, g::Int64 )
    # compute and store no. of individuals in each set
    tracker.ens_pop_dist[g,:] = sum( tracker.pop.sets.h_bar, dims=1 )
    if tracker.verbose_track println("\t $(tracker.ens_pop_dist[g,:])") end
end

function track_running_strat_dist!( tracker::Tracker, g::Int64 )
    # compute and store running average
    if g == 1
        # at g = 1, running average = dist at g = 1 itself
        tracker.ens_running_strat_dist[g,:] = tracker.ens_strat_dist[g,:]
    else
        # at g > 1, running average = (dist at g-1 + dist at g) / g
        tracker.ens_running_strat_dist[g,:] = ( tracker.ens_running_strat_dist[g-1,:]*(g-1) + tracker.ens_strat_dist[g,:] ) / g
    end
    if tracker.verbose_track println("\t $(tracker.ens_running_strat_dist[g,:])") end 
end

function track_opn_dist!( tracker::Tracker, g::Int64 )
    # track opinion distribution
    # opn_base10 = [vectobase3(tracker.pop.sets.h[i,:]) for i in 1:tracker.pop.sets.N] # tracked in sets.jl
    tracker.ens_opn_dist[g,:] = [count(x->x==opn, tracker.pop.sets.h_base10) / tracker.pop.sets.N for opn in tracker.pop.all_opinions] 
    if tracker.verbose_track println("\t $(tracker.ens_opn_dist[g,:])") end 
end

function track_total_means_and_vars!( tracker::Tracker, g::Int64 )
    # track mean and variance of total opinions
    tracker.ens_total_means[g] = mean( sum(tracker.pop.sets.h, dims=2) )
    tracker.ens_total_vars[g]  = var( sum(tracker.pop.sets.h, dims=2) )
    if tracker.verbose_track println("\t $(tracker.ens_total_means[g])") end 
    if tracker.verbose_track println("\t $(tracker.ens_total_vars[g])") end 
end

function track_total_cooperation!( tracker::Tracker, g::Int64, byparty = true )
    # track effective / total cooperation
    # total_interactions::Int64 = 2 * sum([length(x) for x in tracker.pop.sets.set_pairs]) # moved to setup.jl
    tracker.ens_total_cooperation[g] = sum(tracker.pop.prev_actions) / tracker.pop.prev_interactions
    if tracker.pop.sets.P == 2 && byparty

        # track within-party cooperation
        tracker.pop.sim.cooperation  = tracker.pop.prev_actions[:, 1:tracker.pop.sets.N÷2, 1:tracker.pop.sets.N÷2]
        tracker.pop.sim.cooperation += tracker.pop.prev_actions[:, tracker.pop.sets.N÷2+1:tracker.pop.sets.N, tracker.pop.sets.N÷2+1:tracker.pop.sets.N]
        tracker.ens_party_cooperation[g] = sum(tracker.pop.sim.cooperation) / tracker.pop.prev_interactions # can use a different total of needed later

        # track between-party cooperation
        tracker.ens_enemy_cooperation[g] = (sum(tracker.pop.prev_actions) - sum(tracker.pop.sim.cooperation)) / tracker.pop.prev_interactions
    end

    if tracker.verbose_track println("\t $(tracker.ens_total_cooperation[g])") end 
    if tracker.ens_total_cooperation[g] > 1 @warn("!!! total cooperation > 100% !!! set_changed was $(tracker.pop.sets_changed), strategies_changed was $(tracker.pop.strategies_changed)") end
end

function track_opn_simpson!( tracker::Tracker, g::Int64 )
    # track diversity of opinions via simpson's index
    # simpson's index is the probability that two individuals taken at random 
    # from the population represent the same 'type' (here 'type' = set of opinions)
    # opn_base10::Array{Int64, 1} = [vectobase3(tracker.pop.sets.h[i,:]) for i in 1:tracker.pop.sets.N]
    tracker.ens_opn_simpson[g] = sum( [ (count(x->x==opn, tracker.pop.sets.h_base10) / tracker.pop.sets.N)^2 for opn in tracker.pop.all_opinions] )

    if tracker.verbose_track println("\t $(tracker.ens_opn_simpson[g])") end 
end

function track_sets_simpson!( tracker::Tracker, g::Int64 )
    # track diversity of set memberships via simpson index
    # here each 'type' = set of set memberships
    # list_setmems::Array{Int64, 1}   = 1:maximum(tracker.pop.all_opinions) # moved to struct Pop
    # setmems_base10::Array{Int64, 1}= [vectobase3(tracker.pop.sets.h_bar[i,:]) for i in 1:tracker.pop.sets.N] # moved to struct Sets
    tracker.ens_sets_simpson[g] = sum( [ (count(x->x==setmem, tracker.pop.sets.h_bar_base10) / tracker.pop.sets.N)^2 for setmem in tracker.pop.all_memberships] )
    if tracker.verbose_track println("\t $(tracker.ens_sets_simpson[g])") end 
end

function track_cityblock!( tracker::Tracker, g::Int64 )
    # track average pairwise cityblock distance in opinions (h)
    # from Distance.jl (https://github.com/JuliaStats/Distances.jl): cityblock(x, y) = sum(abs(x - y))
    # pairwise computes a symmetric N-by-N matrix of pairwise distances; diagonals are zero
    # compute sum, divide by 2 to remove double counting, and divide by the number of pairs
    tracker.pop.sim.cityblock_ind = pairwise(Cityblock(), tracker.pop.sets.h, dims=1) #h
    tracker.ens_cityblock_ind[g]  = sum( tracker.pop.sim.cityblock_ind ) / 2 / tracker.pop.sets.n_pairs

    # track average cityblock distance between individuals of the same party 
    # pairwise computes a symmetric N/2-by-N/2 matrix of pairwise distances; diagonals are zero
    # compute sum, divide by 2 to remove double counting (can ignore diagonals because they are zero), 
    # divide by 2 again to take average, and divide by the number of pairs per party
    N = tracker.pop.sets.N # to avoid cumbersome notation
    tracker.pop.sim.cityblock_party  = tracker.pop.sim.cityblock_ind[1:N÷2,1:N÷2]     # party 1 only
    tracker.pop.sim.cityblock_party += tracker.pop.sim.cityblock_ind[N÷2+1:N,N÷2+1:N] # party 2 only
    tracker.ens_cityblock_party[g]   = sum(tracker.pop.sim.cityblock_party) / 2 / 2 / tracker.pop.sets.n_party_pairs

    # similarly, track average cityblock distance between individuals of different parties
    # tracker.pop.sim.cityblock_ind[1:N÷2,N÷2+1:N] is the transpose of (same mat)[N÷2+1:N,1:N÷2] so only need to count once
    # instead of complicated divisions, we only need to take the mean of this sub-matrix
    tracker.pop.sim.cityblock_party  = tracker.pop.sim.cityblock_ind[1:N÷2,N÷2+1:N] # party 1 - party 2
    tracker.ens_cityblock_enemy[g]   = mean( tracker.pop.sim.cityblock_party )

end

function track_hamming!( tracker::Tracker, g::Int64 )
    # track average pairwise hamming distance in set memberships (h_bar)
    # from Distance.jl (https://github.com/JuliaStats/Distances.jl): hamming(k, l) = sum(k .!= l)
    # pairwise computes a symmetric N-by-N matrix of pairwise distances; diagonals are zero
    # compute sum, divide by 2 to remove double counting, and divide by the number of pairs
    tracker.pop.sim.hamming_ind = pairwise(Hamming(), tracker.pop.sets.h_bar, dims=1) # h_bar
    tracker.ens_hamming_ind[g]  = sum( tracker.pop.sim.hamming_ind ) / 2 / tracker.pop.sets.n_pairs

    # track average hamming distance between individuals of the same party 
    # pairwise computes a symmetric N/2-by-N/2 matrix of pairwise distances; diagonals are zero
     # compute sum, divide by 2 to remove double counting (can ignore diagonals because they are zero), 
    # divide by 2 again to take average, and divide by the number of pairs per party
    N = tracker.pop.sets.N # to avoid cumbersome notation
    tracker.pop.sim.hamming_party  = tracker.pop.sim.hamming_ind[1:N÷2,1:N÷2]     # party 1 only
    tracker.pop.sim.hamming_party += tracker.pop.sim.hamming_ind[N÷2+1:N,N÷2+1:N] # party 2 only
    tracker.ens_hamming_party[g]   = sum(tracker.pop.sim.hamming_party) / 2 / 2 / tracker.pop.sets.n_party_pairs

    # similarly, track average hamming distance between individuals of different parties
    # tracker.pop.sim.hamming_ind[1:N÷2,N÷2+1:N] is the transpose of (same mat)[N÷2+1:N,1:N÷2] so only need to count once
    # instead of complicated divisions, we only need to take the mean of this sub-matrix
    tracker.pop.sim.hamming_party  = tracker.pop.sim.hamming_ind[1:N÷2,N÷2+1:N] # party 1 - party 2
    tracker.ens_hamming_enemy[g]   = mean( tracker.pop.sim.hamming_party )

end


function make_summary_table( tracker::Tracker, gen_cutoff::Int64, generations::Int64 )
    # function to process data into summary table
    # ignore generations 1 through gen_cutoff

    min_gen = floor(Int64,gen_cutoff)

    short_strat_dist         = tracker.ens_strat_dist[ min_gen:generations,: ]
    short_pop_dist           = tracker.ens_pop_dist[ min_gen:generations,: ]
    short_opn_dist           = tracker.ens_opn_dist[ min_gen:generations,: ]
    short_total_means        = tracker.ens_total_means[ min_gen:generations ]
    short_total_vars         = tracker.ens_total_vars[ min_gen:generations ]
    short_total_cooperation  = tracker.ens_total_cooperation[ min_gen:generations ]
    short_party_cooperation  = tracker.ens_party_cooperation[ min_gen:generations ]
    short_enemy_cooperation  = tracker.ens_enemy_cooperation[ min_gen:generations ]
    short_opn_simpson        = tracker.ens_opn_simpson[ min_gen:generations ]
    short_sets_simpson       = tracker.ens_sets_simpson[ min_gen:generations ]
    short_cityblock_ind      = tracker.ens_cityblock_ind[ min_gen:generations ]
    short_cityblock_party    = tracker.ens_cityblock_party[ min_gen:generations ]
    short_cityblock_enemy    = tracker.ens_cityblock_enemy[ min_gen:generations ]
    short_hamming_ind        = tracker.ens_hamming_ind[ min_gen:generations ]
    short_hamming_party      = tracker.ens_hamming_party[ min_gen:generations ]
    short_hamming_enemy      = tracker.ens_hamming_enemy[ min_gen:generations ]
    short_running_strat_dist = tracker.ens_running_strat_dist[ generations,: ] # just the last number
   
    # put means into a table
    ens_summary = hcat( generations, 
                        tracker.pop.sets.N, tracker.pop.sets.M, tracker.pop.sets.K, tracker.pop.sets.P, 
                        tracker.pop.game.b, tracker.pop.game.c, tracker.pop.game.β, 
                        tracker.pop.game.u, tracker.pop.game.v, tracker.pop.game.p, tracker.pop.game.ϵ,
                        mean(short_strat_dist, dims=1),        var(short_strat_dist, dims=1),
                        mean(short_pop_dist, dims=1),          var(short_pop_dist, dims=1),
                        mean(short_opn_dist, dims=1),          var(short_opn_dist, dims=1),
                        mean(short_total_means, dims=1),       var(short_total_means, dims=1),
                        mean(short_total_vars, dims=1),        var(short_total_vars, dims=1),
                        mean(short_total_cooperation, dims=1), var(short_total_cooperation, dims=1),
                        mean(short_party_cooperation, dims=1), 
                        mean(short_enemy_cooperation, dims=1),
                        mean(short_opn_simpson, dims=1),       var(short_opn_simpson, dims=1),
                        mean(short_sets_simpson, dims=1),      var(short_sets_simpson, dims=1),
                        mean(short_cityblock_ind, dims=1),     var(short_cityblock_ind, dims=1),
                        mean(short_cityblock_party, dims=1),   var(short_cityblock_party, dims=1),
                        mean(short_cityblock_enemy, dims=1),   var(short_cityblock_enemy, dims=1),
                        mean(short_cityblock_ind - short_cityblock_party, dims=1),
                        mean(short_hamming_ind, dims=1),       var(short_hamming_ind, dims=1),
                        mean(short_hamming_party, dims=1),     var(short_hamming_party, dims=1),
                        mean(short_hamming_enemy, dims=1),     var(short_hamming_enemy, dims=1),
                        mean(short_hamming_ind - short_hamming_party, dims=1),
                        short_running_strat_dist'
                        )

    # table_summary = DataFrame(vcat(total_summary...)) # remove later
    col_list = hcat("generations", "N", "M", "K", "P", "b", "c", "β", "u", "v", "p", "ϵ",
                "DD_mean", "DC_mean", "CD_mean", "CC_mean", "DD_var", "DC_var", "CD_var", "CC_var",
                ["set_$(x)_mean" for dummy=1:1, x=1:tracker.pop.sets.M], ["set_$(x)_var" for dummy=1:1, x=1:tracker.pop.sets.M],
                ["opinion_$(x)_mean" for dummy=1:1, x=tracker.pop.all_opinions], 
                ["opinion_$(x)_var" for dummy=1:1, x=tracker.pop.all_opinions], 
                "total_means_mean", "total_means_var", 
                "total_vars_mean", "total_vars_var", 
                "total_cooperation_mean", "total_cooperation_var",
                "party_cooperation_mean", "enemy_cooperation_var",
                "opn_simpson_mean", "opn_simpson_var", 
                "sets_simpson_mean", "sets_simpson_var",
                "cityblock_ind_mean", "cityblock_ind_var", 
                "cityblock_party_mean", "cityblock_party_var", 
                "cityblock_enemy_mean", "cityblock_enemy_var", 
                "cityblock_ind-party_mean",
                "hamming_ind_mean", "hamming_ind_var", 
                "hamming_party_mean", "hamming_party_var", 
                "hamming_enemy_mean", "hamming_enemy_var", 
                "hamming_ind-party_mean", 
                "DD_mean_final", "DC_mean_final", "CD_mean_final", "CC_mean_final"
                )
    colnames = join(col_list, ',')

    return ens_summary, colnames
end

# alternative method
function make_summary_table( tracker::Tracker, gen_cutoff::Int64, generations::Int64, col_list::Array{String, 2} )
    # function to process data into summary table
    # ignore generations 1 through gen_cutoff

    min_gen = floor(Int64,gen_cutoff)

    short_strat_dist         = tracker.ens_strat_dist[ min_gen:generations,: ]
    short_pop_dist           = tracker.ens_pop_dist[ min_gen:generations,: ]
    short_opn_dist           = tracker.ens_opn_dist[ min_gen:generations,: ]
    short_total_means        = tracker.ens_total_means[ min_gen:generations ]
    short_total_vars         = tracker.ens_total_vars[ min_gen:generations ]
    short_total_cooperation  = tracker.ens_total_cooperation[ min_gen:generations ]
    short_party_cooperation  = tracker.ens_party_cooperation[ min_gen:generations ]
    short_enemy_cooperation  = tracker.ens_enemy_cooperation[ min_gen:generations ]
    short_opn_simpson        = tracker.ens_opn_simpson[ min_gen:generations ]
    short_sets_simpson       = tracker.ens_sets_simpson[ min_gen:generations ]
    short_cityblock_ind      = tracker.ens_cityblock_ind[ min_gen:generations ]
    short_cityblock_party    = tracker.ens_cityblock_party[ min_gen:generations ]
    short_cityblock_enemy    = tracker.ens_cityblock_enemy[ min_gen:generations ]
    short_hamming_ind        = tracker.ens_hamming_ind[ min_gen:generations ]
    short_hamming_party      = tracker.ens_hamming_party[ min_gen:generations ]
    short_hamming_enemy      = tracker.ens_hamming_enemy[ min_gen:generations ]
    short_running_strat_dist = tracker.ens_running_strat_dist[ generations,: ] # just the last number
   
    # put means into a table
    ens_summary = hcat( generations, 
                        tracker.pop.sets.N, tracker.pop.sets.M, tracker.pop.sets.K, tracker.pop.sets.P, 
                        tracker.pop.game.b, tracker.pop.game.c, tracker.pop.game.β, 
                        tracker.pop.game.u, tracker.pop.game.v, tracker.pop.game.p, tracker.pop.game.ϵ,
                        mean(short_strat_dist, dims=1),        var(short_strat_dist, dims=1),
                        mean(short_pop_dist, dims=1),          var(short_pop_dist, dims=1),
                        mean(short_opn_dist, dims=1),          var(short_opn_dist, dims=1),
                        mean(short_total_means, dims=1),       var(short_total_means, dims=1),
                        mean(short_total_vars, dims=1),        var(short_total_vars, dims=1),
                        mean(short_total_cooperation, dims=1), var(short_total_cooperation, dims=1),
                        mean(short_party_cooperation, dims=1), 
                        mean(short_enemy_cooperation, dims=1),
                        mean(short_opn_simpson, dims=1),       var(short_opn_simpson, dims=1),
                        mean(short_sets_simpson, dims=1),      var(short_sets_simpson, dims=1),
                        mean(short_cityblock_ind, dims=1),     var(short_cityblock_ind, dims=1),
                        mean(short_cityblock_party, dims=1),   var(short_cityblock_party, dims=1),
                        mean(short_cityblock_enemy, dims=1),   var(short_cityblock_enemy, dims=1),
                        mean(short_cityblock_ind - short_cityblock_party, dims=1),
                        mean(short_hamming_ind, dims=1),       var(short_hamming_ind, dims=1),
                        mean(short_hamming_party, dims=1),     var(short_hamming_party, dims=1),
                        mean(short_hamming_enemy, dims=1),     var(short_hamming_enemy, dims=1),
                        mean(short_hamming_ind - short_hamming_party, dims=1),
                        short_running_strat_dist'
                        )

    colnames = join(col_list, ',')

    return ens_summary, colnames
end

function make_col_list( M::Int64 )
    # function to define the column names for data output

    all_opinions = (-3^M+1)÷2 : (3^M-1)÷2 # list of opinions

    col_list = hcat("generations", "N", "M", "K", "P", "b", "c", "β", "u", "v", "q", "ϵ",
                    "DD_mean", "DC_mean", "CD_mean", "CC_mean", "DD_var", "DC_var", "CD_var", "CC_var",
                    ["set_$(x)_mean" for dummy=1:1, x=1:M], ["set_$(x)_var" for dummy=1:1, x=1:M],
                    ["opinion_$(x)_mean" for dummy=1:1, x=all_opinions], 
                    ["opinion_$(x)_var" for dummy=1:1, x=all_opinions], 
                    "total_means_mean", "total_means_var", 
                    "total_vars_mean", "total_vars_var", 
                    "total_cooperation_mean", "total_cooperation_var",
                    "party_cooperation_mean", "enemy_cooperation_var",
                    "opn_simpson_mean", "opn_simpson_var", 
                    "sets_simpson_mean", "sets_simpson_var",
                    "cityblock_ind_mean", "cityblock_ind_var", 
                    "cityblock_party_mean", "cityblock_party_var", 
                    "cityblock_enemy_mean", "cityblock_enemy_var", 
                    "cityblock_ind-party_mean",
                    "hamming_ind_mean", "hamming_ind_var", 
                    "hamming_party_mean", "hamming_party_var", 
                    "hamming_enemy_mean", "hamming_enemy_var", 
                    "hamming_ind-party_mean", 
                    "DD_mean_final", "DC_mean_final", "CD_mean_final", "CC_mean_final"
                    )

    return col_list
end
