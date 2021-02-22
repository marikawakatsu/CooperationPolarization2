function evolve_neutral!( pop::Population, generations::Int64 = 1 )
    # function to evolve the population under neutral selection
    # that is, role models are selected uniformly at random; fitnesses have no impact
    # generations specifies the total number of generations for which to run the simulation
    for i::Int64 in 1:generations
        if pop.verbose println("\ninitiating generation $(pop.generation)") end
        
        # first, update payoffs
        if pop.verbose println("updating actions only") end
        update_actions_neutral!(pop)  # method that updates actions only
        
        # then, select a learner and a role model and perform update
        if pop.verbose println("") end
        if pop.verbose println("evolving, generation $(pop.generation)") end
        update_strategies_and_opinions_neutral!(pop)
        
        if pop.verbose println("") end
        if pop.verbose println("ending generation $(pop.generation)") end
        pop.generation += 1
    end
end

function update_actions_neutral!( pop::Population )
    # update individual actions
    # for neutral selection (β = 0), so payoffs are not needed
    
    if pop.sets_changed 
        # if (previous) learner's sets have changed, re-compute all actions
        if pop.verbose println("updating all actions because set_changed = $(pop.sets_changed)") end
        pop.prev_actions .= zero(Int64)   # reset all prev_actions
        
        # then, recompute all payoffs
        for k in 1:pop.sets.M  # for each set
            for (i, j) in pop.sets.set_pairs[k]   # for each pair within each set
                compute_actions_neutral!( pop, i, j, k )  # update action between i & j when intereacting in set k
            end
        end
    elseif pop.strategies_changed  
        # if (previous) learner's strategy has changed, recompute actions involving that individual
        if pop.verbose println("updating actions involving $(pop.prev_learner) b/c strategies_changed = $(pop.strategies_changed)") end
        pop.prev_actions[:,pop.prev_learner,:] .= zero(Int64)
        pop.prev_actions[:,:,pop.prev_learner] .= zero(Int64)

        # then, recompute payoffs involving prev_learner
        for k in 1:pop.sets.M  # for each set
            for (i, j) in pop.sets.set_pairs[k]  # for each pair within each set
                if i == pop.prev_learner || j == pop.prev_learner  
                    compute_actions_neutral!( pop, i, j, k )  # update action between i & j when intereacting in set k
                end
            end
        end  
    else 
        if pop.verbose println("not updating actions b/c sets_changed = $(pop.sets_changed) and strategies_changed = $(pop.strategies_changed)") end
        return # no need to update prev_interactions or payoffs
    end
    # sum rows to get total individual payoffs
    pop.prev_interactions = 2 * sum([length(x) for x in pop.sets.set_pairs])

    if pop.verbose println( "\npayoff matrix is now $(pop.payoffs_mat)" ) end
    if pop.verbose println( "payoffs are now $(pop.payoffs)" ) end
    if pop.verbose println( "there were $(pop.prev_interactions) interactions, $(sum(pop.prev_actions)) of which were C" ) end
end

function compute_actions_neutral!( pop::Population, i::Int64, j::Int64, k::Int64 )
    # function to compute payoffs for individuals i and j when they interact in set k
    pop.sim.i_action, pop.sim.j_action = -1, -1  # set to -1 to catch any errors in action assignment

    # actions depend on whether i and j agree on issue k or not
    if pop.sets.h[i,k] == 0 || pop.sets.h[j,k] == 0
        # if either i or j does not belong to set k, something's wrong
        # leave pop.sim.i_action, pop.sim.j_action as -1 so that the compute_payoffs! throws an error when assigning payoffs
        @warn("$i and $j should not be intercting in set $k"); return
    elseif pop.sets.h[i,k] == pop.sets.h[j,k]                              
        pop.sim.i_action, pop.sim.j_action = pop.strategies[i,1], pop.strategies[j,1] # if they agree, play s_agree
    else                                                              
        pop.sim.i_action, pop.sim.j_action = pop.strategies[i,2], pop.strategies[j,2] # if they disagree, play s_disagree
    end

    if pop.verbose
        println("\n\tupdating actions of $i and $j in set $k; agree on issue $k? $(pop.sets.h[i,k] == pop.sets.h[j,k])")
        println("\t\t$i has strats $(pop.strategies[i,:]); $i plays strat $(pop.sim.i_action)")
        println("\t\t$j has strats $(pop.strategies[j,:]); $j plays strat $(pop.sim.j_action)")
    end
    # track actions; no need to track payoffs under neutral selection
    pop.prev_actions[k,i,j] = pop.sim.i_action
    pop.prev_actions[k,j,i] = pop.sim.j_action
end

function update_strategies_and_opinions_neutral!( pop::Population )
    # function to perform death-birth update
    if pop.verbose println("Strategies before updating:\t $(pop.strategies)") end
    if pop.verbose println("Set memberships before updating: $(pop.sets.h)") end
    
    # choose a learner uniformly at random
    pop.sim.learner = sample( 1:pop.sets.N )
    # choose a role model also uniformly at random (for neutral selection)
    pop.sim.role_model = sample( 1:pop.sets.N )

    # learner updates
    if pop.verbose println("\n\trandomly chosen learner is $(pop.sim.learner), role_model is $(pop.sim.role_model)") end
    if pop.verbose println("\t\tall payoffs are $(pop.payoffs)") end
    if pop.verbose println("\t\tfitnesses are $(pop.sim.fitnesses)") end
    
    # compute probability of imitation: 1 if same party, 1-p if different parties
    pop.sets.affils[pop.sim.learner] == pop.sets.affils[pop.sim.role_model] ? imit_prob = 1.0 : imit_prob = 1-pop.game.ps[pop.sets.affils[pop.sim.learner]]
    
    if pop.verbose 
        println("\t\tlearner and role model in the same party? $(pop.sets.affils[ pop.sim.learner ] == pop.sets.affils[ pop.sim.role_model ]), imitate with prob $(imit_prob)") 
    end

    # update strategies and sets: imitation / mutation occurs with probability imit_prob
    if rand() < imit_prob
        if pop.verbose println("\twill consider imitating") end
        update_strategies!( pop ) # includes mutation
        update_opinions!( pop )   # includes mutation
        pop.num_imitating += 1
    else
        if pop.verbose println("\tnah won't imitate") end
        no_strategy_or_opinion_update!( pop )
        pop.num_not_imitating += 1
    end

    # track updates
    pop.prev_learner = pop.sim.learner
    if pop.sets_changed == true
        # h is updated above; need to update h_bar, h_base10, h_bar_base10, set_members, set_pairs
        pop.sets.set_members, pop.sets.set_pairs = set_members_and_pairs(pop.sets.h)  # update set_members and set_pairs
        pop.sets.h_bar        = [abs(x) for x in pop.sets.h]
        pop.sets.h_base10     = [vectobase3(pop.sets.h[i,:]) for i in 1:pop.sets.N] 
        pop.sets.h_bar_base10 = [vectobase3(pop.sets.h_bar[i,:]) for i in 1:pop.sets.N] 
        pop.num_sets_changed += 1
    end
    if pop.strategies_changed == true
        # strategies are updated above; need to update strategies_base10
        pop.strategies_base10 = [2 * pop.strategies[i,1] + pop.strategies[i,2] for i in 1:pop.sets.N]
        pop.num_strategies_changed += 1
    end
    if pop.strategies_changed && pop.sets_changed
        pop.num_both_changed += 1
    end

    if pop.verbose println("\nStrategies after updating:\t $(pop.strategies)") end
    if pop.verbose println("Set memberships after updating:\t $(pop.sets.h)") end
    if pop.verbose println("tracking: prev_learner = $(pop.prev_learner), strategies_changed = $(pop.strategies_changed), sets_changed = $(pop.sets_changed)") end  
end
