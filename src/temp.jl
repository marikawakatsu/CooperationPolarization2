
# scratchwork!
function initialize_interaction_mat( h::Array{Int64, 2} )
    # return a matrix of interactions indicating who interacts with whom
    N, M = size(h)
    interaction_mat = zeros(Int64, N, N, M) # initialize
    for k in 1:M
        for i in 1:N
            if h[i,k] in [-1,1] # check if i belongs to k
                for j in i+1:N
                    if h[j,k] in [-1,1] # check if j both belong to set k
                        interaction_mat[i,j,k] = 1
                        interaction_mat[j,i,k] = 1
                    end
                end
            end
        end
    end
    return interaction_mat
end


function update_interaction_mat!( pop::Population )
    # function to upddate interaction_mat after each round
    if pop.sets_changed
        pop.sets.interaction_mat .= 0
        for k = 1:M
            for i = 1:N
                if h[i,k] in [-1,1] # check if i belongs to k
                    for j = i+1:N
                        if h[j,k] in [-1,1] # check if j both belong to set k
                            pop.sets.interaction_mat[i,j,k] = 1
                            pop.sets.interaction_mat[j,i,k] = 1
                        end
                    end
                end
            end
        end
    end
end

function update_interaction_mat2!( pop::Population )
    # function to upddate interaction_mat after each round
    if pop.sets_changed
        # only update interaction matrix if a set change has occurred
        set_members = [Int64[] for k in 1:M]
        for indices in findall(x->x in [1,-1], pop.sets.h) 
            push!(set_members[indices[2]],indices[1]) 
        end
        
        # update interaction matrix
        pop.sets.interaction_mat .= 0
        for k = 1:M
            for i in set_members[k], j in set_members[k]
                if i < j
                    pop.sets.interaction_mat[i,j,k] = 1
                    pop.sets.interaction_mat[j,i,k] = 1
                end
            end
        end
    end
end

for k in 1:pop.sets.M  # for each set
    for (i, j) in pop.sets.set_pairs[k]   # for each pair within each set
        compute_payoffs!( pop, i, j, k )  # compute payoffs for i & j when intereacting in set k
                                          # this also updates pop.prev_actions accordingly
    end
end

for k in 1:pop.sets.M  
    for i in 1:pop.sets.N
        for indices in findall(x->x==1, pop.sets.interaction_mat[i,:,k])
            # fill in here
        end
    end
end

