# simulate.jl: functions to simulate the model
# include("setup.jl")
function evolve!( pop::Population, generations::Int64 = 1 )
    # from ReputationSets https://github.com/tkessinger/reputation_sets
    # the main function to evolve the population
    # generations specifies the total number of generations for which to run the simulation
    for i::Int64 in 1:generations
        if pop.verbose println("\ninitiating generation $(pop.generation)") end
        
        # first, update payoffs
        if pop.verbose println("updating actions and payoffs") end
        update_payoffs_using_matrix!(pop)  # method that stores payoffs in a matrix
        
        # then, select a learner and a role model and perform update
        if pop.verbose println("") end
        if pop.verbose println("evolving, generation $(pop.generation)") end
        update_strategies_and_opinions_db!(pop)
        
        if pop.verbose println("") end
        if pop.verbose println("ending generation $(pop.generation)") end
        pop.generation += 1
    end
end

function update_payoffs_using_matrix!( pop::Population )
    # update individual payoffs by computing pairwise payoffs and storing them in a matrix
    
    # to make the code run faster, we only update payoffs in a given generation 
    # if strategies or sets have changed. sets_changed and strategies_changed are 
    # initialized as true so that payoffs are computed in the first round

    if pop.sets_changed 
        # if (previous) learner's sets have changed, re-compute all payoffs
        if pop.verbose println("updating all payoffs because set_changed = $(pop.sets_changed)") end
        pop.payoffs_mat  .= zero(Float64) # reset all payoffs_mat
        pop.prev_actions .= zero(Float64) # reset all prev_actions
        
        # then, recompute all payoffs
        for k in 1:pop.sets.M  # for each set
            for (i, j) in pop.sets.set_pairs[k]   # for each pair within each set
                compute_payoffs!( pop, i, j, k )  # compute payoffs for i & j when intereacting in set k
                                                  # this also updates pop.prev_actions accordingly
            end
        end
    elseif pop.strategies_changed  
        # if (previous) learner's strategy has changed, recompute payoffs involving that individual
        if pop.verbose println("updating payoffs involving $(pop.prev_learner) because strategies_changed = $(pop.strategies_changed)") end
        pop.payoffs_mat[pop.prev_learner,:]    .= zero(Float64) # first, reset payoffs involving prev_learner
        pop.payoffs_mat[:,pop.prev_learner]    .= zero(Float64)
        pop.prev_actions[:,pop.prev_learner,:] .= zero(Float64)
        pop.prev_actions[:,:,pop.prev_learner] .= zero(Float64)

        # then, recompute payoffs involving prev_learner
        for k in 1:pop.sets.M  # for each set
            for (i, j) in pop.sets.set_pairs[k]  # for each pair within each set
                if i == pop.prev_learner || j == pop.prev_learner  
                    compute_payoffs!( pop, i, j, k )  # compute payoffs for i & j when intereacting in set k
                                                      # this also updates pop.prev_actions accordingly
                end
            end
        end      
    else 
        if pop.verbose println("not updating payoffs b/c sets_changed = $(pop.sets_changed) and strategies_changed = $(pop.strategies_changed)") end
        return # no need to update prev_interactions or payoffs
    end

    # sum rows to get total individual payoffs
    pop.payoffs = [round( x, digits=13 ) for x in vec( sum(pop.payoffs_mat, dims=2) )]  # rounding avoids precision errors
    pop.prev_interactions = 2 * sum([length(x) for x in pop.sets.set_pairs])

    if pop.verbose println( "payoff matrix is now $(pop.payoffs_mat)" ) end
    if pop.verbose println( "payoffs are now $(pop.payoffs)" ) end
    if pop.verbose println( "there were $(pop.prev_interactions) interactions, $(sum(pop.prev_actions)) of which were C" ) end
end

function compute_payoffs!( pop::Population, i::Int64, j::Int64, k::Int64 )
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
        println("\tupdating actions of $i and $j in set $k; agree on issue $k? $(pop.sets.h[i,k] == pop.sets.h[j,k])")
        println("\t\t$i has strats $(pop.strategies[i,:]); $i plays strat $(pop.sim.i_action), earns a payoff of $(pop.game.A[pop.sim.i_action+1, pop.sim.j_action+1])")
        println("\t\t$j has strats $(pop.strategies[j,:]); $j plays strat $(pop.sim.j_action), earns a payoff of $(pop.game.A[pop.sim.j_action+1, pop.sim.i_action+1])")
    end

    # add i and j's payoffs to the payoff matrix according to their strategies and the game matrix
    pop.payoffs_mat[i,j] += pop.game.A[pop.sim.i_action+1, pop.sim.j_action+1] # +1 because julia is 1-indexed
    pop.payoffs_mat[j,i] += pop.game.A[pop.sim.j_action+1, pop.sim.i_action+1]

    # track actions
    pop.prev_actions[k,i,j] = pop.sim.i_action
    pop.prev_actions[k,j,i] = pop.sim.j_action
end

function update_strategies_and_opinions_db!( pop::Population )
    # function to perform death-birth update
    if pop.verbose println("Strategies before updating:\t $(pop.strategies)") end
    if pop.verbose println("Set memberships before updating: $(pop.sets.h)") end
    
    # randomly choose a learner and store a list of non_learners
    pop.sim.learner      = sample(1:pop.sets.N)
    pop.sim.non_learners = filter(x->x!=pop.sim.learner, 1:pop.sets.N)
    
    # compute the fitnesses of every other individual in the population
    # pop.sim.fitnesses = 1.0 .+ pop.game.β * pop.payoffs[ pop.sim.non_learners ]
    pop.sim.fitnesses = [1.0 + pop.game.β * pop.payoffs[i] for i in pop.sim.non_learners ] # faster
    if any(x->x<0.0, pop.sim.fitnesses) @error("!!! negative fitness detected, aborting !!!") end
    
    # choose a role model from the rest of the population, weighted by fitness
    pop.sim.role_model = sample( pop.sim.non_learners, Weights(pop.sim.fitnesses) )
    
    # learner updates
    if pop.verbose println("\n\trandomly chosen learner is $(pop.sim.learner), role_model is $(pop.sim.role_model)") end
    if pop.verbose println("\t\tall payoffs are $(pop.payoffs)") end
    if pop.verbose println("\t\tfitnesses (learner excluded) are $(pop.sim.fitnesses)") end
    
    # compute probability of imitation: 1 if same party, q if different parties
    pop.sets.affiliations[pop.sim.learner] == pop.sets.affiliations[pop.sim.role_model] ? imit_prob = 1.0 : imit_prob = pop.game.q
    if pop.verbose 
        println("\t\tlearner and role model in the same party? $(pop.sets.affiliations[ pop.sim.learner ] == pop.sets.affiliations[ pop.sim.role_model ])") 
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
    if pop.verbose 
        println("tracking: prev_learner = $(pop.prev_learner), strategies_changed = $(pop.strategies_changed), sets_changed = $(pop.sets_changed)") 
    end  
end

function update_strategies!( pop::Population )
    # function to update strategy of learner based on role_model
    pop.sim.strategy_copy = copy(pop.strategies[pop.sim.learner,:])  # save a copy to track changes
    # copy suffices here because the assignments below replace elements of pop.sets.h with Int64 (immutable)

    if rand() < pop.game.u
        # strategy mutation occurs
        pop.strategies[pop.sim.learner,:] = rand(collect(0:1), 1, 2) # conditional case
        pop.num_random_strategies += 1
        if pop.verbose println("\t\tmutating $(pop.sim.learner) to strategy $(pop.strategies[pop.sim.learner,:])") end
    else
        # learner imitates the strategy of role model
        pop.strategies[pop.sim.learner,:] = pop.strategies[pop.sim.role_model,:]
        if pop.verbose println("\t\tadopting $(pop.sim.role_model)'s strategy $(pop.strategies[pop.sim.role_model,:])") end
    end
    pop.strategies_changed = (pop.sim.strategy_copy != pop.strategies[pop.sim.learner,:]) # track changes
end

function update_opinions!( pop::Population )
    # function to update set memberships / opinions of learner based on role_model
    pop.sim.sets_copy = copy(pop.sets.h[pop.sim.learner,:])     # save a copy to track changes
    # copy suffices here because the assignments below replace elements of pop.sets.h with Int64 (immutable)
    
    if rand() < pop.game.v 
        # with probability v, set mutation occurs
        if rand() < pop.game.gamma * pop.game.q
            # with probability q*gamma
            # select a random arrangement of K opinions
            pop.sets.h[pop.sim.learner,:] = shuffle( vcat( zeros(Int64, pop.sets.M-pop.sets.K), rand([-1,1], pop.sets.K) ) )
            pop.num_random_opinions += 1
            if pop.verbose println("\t\tmutating $(pop.sim.learner) to random opinions $(pop.sets.h[pop.sim.learner,:])") end
        else
            # with probability 1-q*gamma
            # select a random set of K issues, select opinions corresponding to party
            if pop.sets.affiliations[pop.sim.learner] == 1
                pop.sets.h[pop.sim.learner,:] = shuffle( vcat( zeros(Int64, pop.sets.M-pop.sets.K), -ones(Int64, pop.sets.K) ) )
            elseif pop.sets.affiliations[pop.sim.learner] == 2
                pop.sets.h[pop.sim.learner,:] = shuffle( vcat( zeros(Int64, pop.sets.M-pop.sets.K), ones(Int64, pop.sets.K) ) )
            end
            pop.num_biased_opinions += 1
            if pop.verbose 
                println("\t\t$(pop.sim.learner) is in party $(pop.sets.affiliations[ pop.sim.learner ]), mutating to biased opinions $(pop.sets.h[ pop.sim.learner,: ])") 
            end
        end
    else
        # with probability 1-v, learner imitates the strategy of role model
        pop.sets.h[ pop.sim.learner,: ] = pop.sets.h[ pop.sim.role_model,: ]
        if pop.verbose println("\t\tadopting $( pop.sim.role_model )'s set membership $(pop.sets.h[ pop.sim.role_model,: ])") end
    end
    # h is updated above; need to update h_bar, sets_changed, set_members, set_pairs
    # pop.sets.set_members, pop.sets.set_pairs = set_members_and_pairs(pop.sets.h)  # update set_members and set_pairs # moved up
    pop.sets_changed = (pop.sim.sets_copy != pop.sets.h[ pop.sim.learner,: ]) # track changes
end

function no_strategy_or_opinion_update!( pop::Population )
    # function to keep track of changes when no imitation/mutation occurs
    pop.strategies_changed = false  # track changes
    pop.sets_changed       = false
end