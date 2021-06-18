mutable struct NeutralTracker
    # structure for storing simulation data
    pop::Population
    generations::Int64
    verbose_track::Bool                         # print messages as data are tracked
    ens_strat_dist::Array{Float64, 2}           # track frequencies of strategies
    ens_strat_running::Array{Float64, 2}        # track running average of frequencies of strategies

    ens_y::Array{Float64, 1}                    # track quantity y = P(s_ia == s_ja | i != j)
    ens_y_running::Array{Float64, 1}            # track running average of y
    ens_z::Array{Float64, 1}                    # track quantity z = <h_i h_l | i != j>
    ens_z_running::Array{Float64, 1}            # track running average of z
    ens_g::Array{Float64, 1}                    # track quantity g = <h_i h_l I(s_ia == s_la) | i != l>
    ens_g_running::Array{Float64, 1}            # track running average of g
    ens_h::Array{Float64, 1}                    # track quantity h = <h_j h_l I(s_ia == s_ja) | i != j != l != i>
    ens_h_running::Array{Float64, 1}            # track running average of h

    ens_zp::Array{Float64, 1}                    # track quantity z' = < |h_i||h_l| | i != j>
    ens_zp_running::Array{Float64, 1}            # track running average of z'
    ens_gp::Array{Float64, 1}                    # track quantity g' = < |h_i||h_l| I(s_ia == s_la) | i != l>
    ens_gp_running::Array{Float64, 1}            # track running average of g'
    ens_hp::Array{Float64, 1}                    # track quantity h' = < |h_i||h_l| I(s_ia == s_ja) | i != j != l != i>
    ens_hp_running::Array{Float64, 1}            # track running average of h'

    ens_sia_sjd::Array{Float64, 1}              # track quantity q = <s_ia s_jd | i, j> # not requiring i != j
    ens_sia_sjd_running::Array{Float64, 1}      # track running average of <s_ia s_jd | i, j> # not requiring i != j

    # other options
    ens_sia_sid::Array{Float64, 1}              # track quantity <s_ia s_id> 
    ens_sia_sid_running::Array{Float64, 1}      # track quantity <s_ia s_id> 
    # ens_bracket_sia_sjd::Array{Float64, 1}   # track quantity <s_ia s_jd>
    # ens_bracket_sia_sja::Array{Float64, 1}   # track quantity <s_ia s_ja>
    # ens_bracket_sid_sjd::Array{Float64, 1}   # track quantity <s_id s_jd> # should be identical to above
    # ens_bracket_sia_hi_hl::Array{Float64, 1}   # track quantity <s_ia h_i h_l>
    # ens_bracket_sia_sid_hi_hl::Array{Float64, 1}   # track quantity <s_ia s_la h_i h_l>
    # ens_bracket_sia_sla_hi_hl::Array{Float64, 1}   # track quantity <s_ia s_la h_i h_l>
    # ens_bracket_sia_sld_hi_hl::Array{Float64, 1}   # track quantity <s_ia s_ld h_i h_l>
    # ens_bracket_sia_sja_hj_hl::Array{Float64, 1}   # track quantity <s_ia s_ja h_i h_l> # triplet
    # ens_bracket_sia_sja_hj_hl::Array{Float64, 1}   # track quantity <s_ia s_ja h_i h_l> # triplet

    # constructor
    function NeutralTracker( pop::Population, generations::Int64, verbose_track::Bool )
        # initialize n_trackers
        ens_strat_dist    = Float64[0.0 for gens=1:generations, strats=1:length(pop.all_strategies)]
        ens_strat_running = Float64[0.0 for gens=1:generations, strats=1:length(pop.all_strategies)]

        ens_y         = Float64[0.0 for gens=1:generations]
        ens_y_running = Float64[0.0 for gens=1:generations]
        ens_z         = Float64[0.0 for gens=1:generations]
        ens_z_running = Float64[0.0 for gens=1:generations]
        ens_g         = Float64[0.0 for gens=1:generations]
        ens_g_running = Float64[0.0 for gens=1:generations]
        ens_h         = Float64[0.0 for gens=1:generations]
        ens_h_running = Float64[0.0 for gens=1:generations]

        ens_zp         = Float64[0.0 for gens=1:generations]
        ens_zp_running = Float64[0.0 for gens=1:generations]
        ens_gp         = Float64[0.0 for gens=1:generations]
        ens_gp_running = Float64[0.0 for gens=1:generations]
        ens_hp         = Float64[0.0 for gens=1:generations]
        ens_hp_running = Float64[0.0 for gens=1:generations]

        ens_sia_sjd         = Float64[0.0 for gens=1:generations]
        ens_sia_sjd_running = Float64[0.0 for gens=1:generations]
        ens_sia_sid         = Float64[0.0 for gens=1:generations]
        ens_sia_sid_running = Float64[0.0 for gens=1:generations]

        return new(pop, generations, verbose_track,
                ens_strat_dist, ens_strat_running,
                ens_y, ens_y_running, ens_z, ens_z_running,
                ens_g, ens_g_running, ens_h, ens_h_running,
                ens_zp, ens_zp_running, ens_gp, ens_gp_running, 
                ens_hp, ens_hp_running,
                ens_sia_sjd, ens_sia_sjd_running, 
                ens_sia_sid, ens_sia_sid_running
                )
    end
end

function track_neutral!( n_tracker::NeutralTracker, gen::Int64 )
    # main function to track the population dynamics under neutral selection

    if n_tracker.verbose_track println("generation $gen") end

    track_strat_dist!(n_tracker, gen) # track strategy distribution and running average 

    # singles
    for i in 1:n_tracker.pop.sets.N n_tracker.ens_sia_sid[gen] += get_sia_times_sid(n_tracker, i) end
    n_tracker.ens_sia_sid[gen] = n_tracker.ens_sia_sid[gen] / n_tracker.pop.sets.N

    # duplets, i != j
    for i in 1:n_tracker.pop.sets.N # keep it simple
        for j in i+1:n_tracker.pop.sets.N # i != j
            n_tracker.ens_y[gen] += get_sia_eq_sja(n_tracker, i, j)
            n_tracker.ens_z[gen] += get_hi_dot_hj(n_tracker, i, j)
            n_tracker.ens_g[gen] += get_hi_dot_hj_sia_eq_sja(n_tracker, i, j)

            n_tracker.ens_zp[gen] += get_hbari_dot_hbarj(n_tracker, i, j)
            n_tracker.ens_gp[gen] += get_hbari_dot_hbarj_sia_eq_sja(n_tracker, i, j)
        end
    end
    n_tracker.ens_y[gen]  = n_tracker.ens_y[gen] / n_tracker.pop.sets.n_pairs
    n_tracker.ens_z[gen]  = n_tracker.ens_z[gen] / n_tracker.pop.sets.n_pairs
    n_tracker.ens_g[gen]  = n_tracker.ens_g[gen] / n_tracker.pop.sets.n_pairs
    n_tracker.ens_zp[gen] = n_tracker.ens_zp[gen] / n_tracker.pop.sets.n_pairs
    n_tracker.ens_gp[gen] = n_tracker.ens_gp[gen] / n_tracker.pop.sets.n_pairs

    # duplets, i == j included
    for i in 1:n_tracker.pop.sets.N # keep it simple
        for j in i:n_tracker.pop.sets.N # i == j in included
            n_tracker.ens_sia_sjd[gen] += get_sia_times_sjd(n_tracker, i, j)
        end
    end
    n_tracker.ens_sia_sjd[gen] = n_tracker.ens_sia_sjd[gen] / (n_tracker.pop.sets.n_pairs + n_tracker.pop.sets.N)

    # triplets, i != j != l (!= i)
    triplet_counter = 0 # eventually store this in the n_tracker struct somewhere
    for i in 1:n_tracker.pop.sets.N # keep it simple
        for j in i+1:n_tracker.pop.sets.N # i != j
            for l in j+1:n_tracker.pop.sets.N
                n_tracker.ens_h[gen] += get_hi_dot_hl_sia_eq_sja(n_tracker, i, j, l)
                n_tracker.ens_h[gen] += get_hi_dot_hl_sia_eq_sja(n_tracker, j, l, i)
                n_tracker.ens_h[gen] += get_hi_dot_hl_sia_eq_sja(n_tracker, l, i, j)

                n_tracker.ens_hp[gen] += get_hbari_dot_hbarl_sia_eq_sja(n_tracker, i, j, l)
                n_tracker.ens_hp[gen] += get_hbari_dot_hbarl_sia_eq_sja(n_tracker, j, l, i)
                n_tracker.ens_hp[gen] += get_hbari_dot_hbarl_sia_eq_sja(n_tracker, l, i, j)

                triplet_counter += 3
            end
        end
    end
    n_tracker.ens_h[gen]  = n_tracker.ens_h[gen] / triplet_counter
    n_tracker.ens_hp[gen] = n_tracker.ens_hp[gen] / triplet_counter

    # update running average
    if gen == 1 
        # at gen = 1, running average = dist at gen = 1 itself
        n_tracker.ens_y_running[gen] = n_tracker.ens_y[gen]
        n_tracker.ens_z_running[gen] = n_tracker.ens_z[gen]
        n_tracker.ens_g_running[gen] = n_tracker.ens_g[gen]
        n_tracker.ens_h_running[gen] = n_tracker.ens_h[gen]
        n_tracker.ens_zp_running[gen] = n_tracker.ens_zp[gen]
        n_tracker.ens_gp_running[gen] = n_tracker.ens_gp[gen]
        n_tracker.ens_hp_running[gen] = n_tracker.ens_hp[gen]
        n_tracker.ens_sia_sjd_running[gen] = n_tracker.ens_sia_sjd[gen]
        n_tracker.ens_sia_sid_running[gen] = n_tracker.ens_sia_sid[gen]
    else 
        # at gen > 1, running average = ( (running avg y at gen-1) * (gen-1) + y at gen) / gen
        n_tracker.ens_y_running[gen] = ( n_tracker.ens_y_running[gen-1]*(gen-1) + n_tracker.ens_y[gen] ) / gen
        n_tracker.ens_z_running[gen] = ( n_tracker.ens_z_running[gen-1]*(gen-1) + n_tracker.ens_z[gen] ) / gen
        n_tracker.ens_g_running[gen] = ( n_tracker.ens_g_running[gen-1]*(gen-1) + n_tracker.ens_g[gen] ) / gen
        n_tracker.ens_h_running[gen] = ( n_tracker.ens_h_running[gen-1]*(gen-1) + n_tracker.ens_h[gen] ) / gen
        n_tracker.ens_zp_running[gen] = ( n_tracker.ens_zp_running[gen-1]*(gen-1) + n_tracker.ens_zp[gen] ) / gen
        n_tracker.ens_gp_running[gen] = ( n_tracker.ens_gp_running[gen-1]*(gen-1) + n_tracker.ens_gp[gen] ) / gen
        n_tracker.ens_hp_running[gen] = ( n_tracker.ens_hp_running[gen-1]*(gen-1) + n_tracker.ens_hp[gen] ) / gen
        n_tracker.ens_sia_sjd_running[gen] = ( n_tracker.ens_sia_sjd_running[gen-1]*(gen-1) + n_tracker.ens_sia_sjd[gen] ) / gen
        n_tracker.ens_sia_sid_running[gen] = ( n_tracker.ens_sia_sid_running[gen-1]*(gen-1) + n_tracker.ens_sia_sid[gen] ) / gen
    end
end

function get_sia_times_sid( n_tracker::NeutralTracker, i::Int64 )
    # for <s_ia s_id>
    return Float64( n_tracker.pop.strategies[i,1] * n_tracker.pop.strategies[i,2] )
end

function get_sia_eq_sja( n_tracker::NeutralTracker, i::Int64, j::Int64 )
    # for y = P(s_ia == s_ja | i != j)
    return Float64( n_tracker.pop.strategies[i,1] == n_tracker.pop.strategies[j,1] )
end

function get_sia_times_sjd( n_tracker::NeutralTracker, i::Int64, j::Int64 )
    # for q = <s_ia s_jd | i, j>, including i == j
    return Float64( n_tracker.pop.strategies[i,1] * n_tracker.pop.strategies[j,2] )
end

function get_hi_dot_hj( n_tracker::NeutralTracker, i::Int64, j::Int64 )
    # for z = <h_i h_l | i != j>
    product = 0.
    for k in 1:n_tracker.pop.sets.M
        product += Float64( n_tracker.pop.sets.h[i,k] * n_tracker.pop.sets.h[j,k] )
    end
    return product
    # return Float64( n_tracker.pop.sets.h[i] * n_tracker.pop.sets.h[j] )
end

function get_hi_dot_hj_sia_eq_sja( n_tracker::NeutralTracker, i::Int64, j::Int64 )
    # for g = <h_i h_l I(s_ia == s_la) | i != l>
    # ie h_i h_l if s_ia = s_la is true, else 0
    product = 0.
    if n_tracker.pop.strategies[i,1]==n_tracker.pop.strategies[j,1]
        for k in 1:n_tracker.pop.sets.M
            product += Float64( n_tracker.pop.sets.h[i,k] * n_tracker.pop.sets.h[j,k] )
        end
    end
    return product
    # return n_tracker.pop.strategies[i,1]==n_tracker.pop.strategies[j,1] ? Float64( n_tracker.pop.sets.h[i] * n_tracker.pop.sets.h[j] ) : 0.0
end

function get_hi_dot_hl_sia_eq_sja( n_tracker::NeutralTracker, i::Int64, j::Int64, l::Int64)
    # for h = <h_j h_l I(s_ia == s_ja) | i != j != l != i>
    # ie h_i h_l if s_ia = s_ja is true, else 0
    product = 0.
    if n_tracker.pop.strategies[i,1]==n_tracker.pop.strategies[j,1]
        for k in 1:n_tracker.pop.sets.M
            product += Float64( n_tracker.pop.sets.h[i,k] * n_tracker.pop.sets.h[l,k] )
        end
    end
    return product
    # return n_tracker.pop.strategies[i,1]==n_tracker.pop.strategies[j,1] ? Float64( n_tracker.pop.sets.h[i] * n_tracker.pop.sets.h[l] ) : 0.0
end

function get_hbari_dot_hbarj( n_tracker::NeutralTracker, i::Int64, j::Int64 )
    # for z' = < |h_i||h_l| | i != j>
    product = 0.
    for k in 1:n_tracker.pop.sets.M
        product += Float64( n_tracker.pop.sets.h_bar[i,k] * n_tracker.pop.sets.h_bar[j,k] )
    end
    return product
end

function get_hbari_dot_hbarj_sia_eq_sja( n_tracker::NeutralTracker, i::Int64, j::Int64 )
    # for g' = < |h_i||h_l| I(s_ia == s_la) | i != l>
    product = 0.
    if n_tracker.pop.strategies[i,1]==n_tracker.pop.strategies[j,1]
        for k in 1:n_tracker.pop.sets.M
            product += Float64( n_tracker.pop.sets.h_bar[i,k] * n_tracker.pop.sets.h_bar[j,k] )
        end
    end
    return product
end

function get_hbari_dot_hbarl_sia_eq_sja( n_tracker::NeutralTracker, i::Int64, j::Int64, l::Int64)
    # for h' = < |h_j||h_l| I(s_ia == s_ja) | i != j != l != i>
    product = 0.
    if n_tracker.pop.strategies[i,1]==n_tracker.pop.strategies[j,1]
        for k in 1:n_tracker.pop.sets.M
            product += Float64( n_tracker.pop.sets.h_bar[i,k] * n_tracker.pop.sets.h_bar[l,k] )
        end
    end
    return product
end

function track_strat_dist!( n_tracker::NeutralTracker, gen::Int64 )
    # compute and store frequency of individuals with each strategy
    # strategies_base10           = [2*tracker.pop.strategies[i,1] + n_tracker.pop.strategies[i,2] for i in 1:tracker.pop.sets.N]
    n_tracker.ens_strat_dist[gen,:] = [count(x->x==strat, n_tracker.pop.strategies_base10) / n_tracker.pop.sets.N for strat in n_tracker.pop.all_strategies] 

    # track running average
    if gen == 1 # at g = 1, running average = dist at g = 1 itself
        n_tracker.ens_strat_running[gen,:] = n_tracker.ens_strat_dist[gen,:]
    else # at g > 1, running average = (dist at g-1 + dist at g) / g
        n_tracker.ens_strat_running[gen,:] = ( n_tracker.ens_strat_running[gen-1,:]*(gen-1) + n_tracker.ens_strat_dist[gen,:] ) / gen
    end
end