# Setup.jl: define structs Sets, Game, and Population
# __precompile__()
# module Setup

export Sets, Game, Population
export compute_total_memberships
export party_assignment
export set_members_and_pairs
export random_sets

# load packages
using Base
using Random
using StatsBase
using Statistics
using Combinatorics

mutable struct Sets
    # structure for storing information about set memberships, opinions, and parties

    N::Int64 # number of individuals
    M::Int64 # number of sets
    K::Int64 # number of sets to which each individual belongs
    P::Int64 # number of parties to which each individual belongs
    set_members::Array{Array{Int64, 1}, 1}            # lists of individuals in each set
    set_pairs::Array{Array{Tuple{Int64, Int64}, 1},1} # list of tuples, each tuple is a pair of individuals in a set
    h::Array{Int64, 2}      # N-by-M int array denoting *opinions*
                            # call h[:, k] to get all opinions in set k
                            # or h[i, :] to get all of i's opinions
    h_bar::Array{Int64,2}   # N-by-M int array denoting *set memberships*
                            # call h_bar[:, k] to get all members in set k
                            # or h_bar[i, :] to get all of i's set memberships
    h_base10::Array{Int64, 1}      # N-by-1 vector denoting *opinions* in base 10 (from -(3^M-1)/2 to (3^M-1)/2)
    h_bar_base10::Array{Int64, 1}  # N-by-1 vector denoting *set memberships* in base 10 (from 1 to (3^M-1)/2)
    affiliations::Array{Int64, 1}  # N-by-1 int array denoting *party affiliations* (0, ..., P-1)
                            # call p[i] to get part membership of i
    n_pairs::Int64          # number of pairs = binomial(N,2) for computing pairwise distances later
    n_party_pairs::Int64    # number of within-party pairs = binomial(N÷2, 2)
    
    # constructors for pre-specified h
    function Sets( h::Array{Int64, 2}, P::Int64=2 )
        
        N, M    = size(h) # extract numbers of individuals and sets
        h_bar   = [abs(x) for x in h]              # extract set memberships from h
        K       = compute_total_memberships(h_bar) # each member must belong to K sets
        affiliations           = party_assignment(N,P)            # each member is affiliated with a party
        n_pairs, n_party_pairs = binomial(N, 2), binomial(N÷2, 2) # number of pairs
        set_members, set_pairs = set_members_and_pairs(h_bar)   

        h_base10     = [vectobase3(h[i,:]) for i in 1:N] # track opinions in base10
        h_bar_base10 = [vectobase3(h_bar[i,:]) for i in 1:N] # track set memberships in base10

        return new( N, M, K, P, set_members, set_pairs, h, h_bar, h_base10, h_bar_base10, affiliations, n_pairs, n_party_pairs)
    end   
    
    # alternative constructor with N, M, K given explicitly
    function Sets( N::Int64, M::Int64, K::Int64, h::Array{Int64, 2}, P::Int64=2 )
            
        h_bar                  = [abs(x) for x in h]              # extract set memberships from h
        affiliations           = party_assignment(N,P)            # each member belongs to a party
        n_pairs, n_party_pairs = binomial(N, 2), binomial(N÷2, 2) # number of pairs
        set_members, set_pairs = set_members_and_pairs(h_bar)  

        h_base10     = [vectobase3(h[i,:]) for i in 1:N] # track opinions in base10
        h_bar_base10 = [vectobase3(h_bar[i,:]) for i in 1:N] # track set memberships in base10

        return new( N, M, K, P, set_members, set_pairs, h, h_bar, h_base10, h_bar_base10, affiliations, n_pairs, n_party_pairs)
    end
end

function compute_total_memberships( h_bar::Array{Int64, 2} )::Int64
    # function to compute K, the number of sets to which each individual belongs
    total_memberships = sum( h_bar, dims=2 ) # total no. of sets indivduals belong to
    # check that all individuals belong to the same number of sets
    if length( unique(total_memberships) ) == 1   
        K = total_memberships[1] # K is the unique element of total_memberships
    else
        @warn("check that each individual belongs to K sets"); return 
    end

    return K
end

function set_members_and_pairs( h::Array{Int64, 2} )
    # return a list of list of set members given a boolean array of memberships

    N, M = size(h)
    set_members = [Int64[] for k in 1:M]
    for indices in findall(x->x in [1,-1], h) 
        push!(set_members[indices[2]],indices[1]) 
    end
    
    set_pairs = [Tuple{Int64, Int64}[] for k in 1:M]
    for k in 1:M
        for i in set_members[k]
            for j in set_members[k]
                if i < j push!(set_pairs[k], (i, j)) end
            end
        end
    end
    return set_members, set_pairs
end

function party_assignment( N::Int64, P::Int64=2 )
    # return party assignments p given N individuals and P parties
    # P = 2 only for now
    if N % P != 0
        @warn("Parties are uneven in size!")
    elseif P == 1
        # everyone is in the same party (party 1)
        affiliations = ones(Int64, N)
    elseif P == 2
        # without loss of generality, first N/2 are in party 1, second N/2 in party 2
        # party 0 is reserved for independents
        affiliations = vcat( ones(Int64, N÷2), repeat([2],N÷2) ) 
    elseif P != 2
        @warn("Party assignment for P>2 is undefined")
    end
    return affiliations
end

function random_sets( N::Int64, M::Int64, K::Int64, P::Int64=2 )
    # cosntructor for randomized set membership
    # where each individual belongs to K > 0 sets (no loners)
    h = zeros(Int64, N, M) # initializing

    if P == 2
        # for each individual, pick K random sets to belong to
        # for each of the K sets, choose -1 for party 1 and +1 for party 2
        for i in 1:N÷2   h[i,:] = shuffle(vcat(zeros(Int64,M-K),-ones(Int64,K))) end # party 1
        for i in N÷2+1:N h[i,:] = shuffle(vcat(zeros(Int64,M-K), ones(Int64,K))) end # party 2
    else
        # for each individual, pick K random sets to belong to
        # for each of the K sets, choose -1 or +1 opinion randomly
        for i in 1:N h[i,:] = shuffle( vcat(zeros(Int64,M-K),rand([-1,1],K)) ) end
    end
    return Sets(h)
end

function vectobase3( vec::Array{Int64, 1} )::Int64
    # convert vector of opinions to int 
    # by interpreting the vector as digits in base 3
    foldl( (x::Int64, y::Int64) -> 3 * x + y, vec )
end

struct Game 
    # structure for storing game parameters

    b::Float64    # benefit to cooperating
    c::Float64    # cost to cooperating
    β::Float64    # selection strength
    u::Float64    # strategy mutation rate
    v::Float64    # set mutation rate
    p::Float64    # party bias (p = 0 no bias, p = 1 strong bias)
    ε::Float64    # damping parameter attenuating party bias (ε = 1 no damping, ε = 0 full damping)
    A::Array{Float64, 2} # the game matrix

    function Game( b::Float64, c::Float64, β::Float64, u::Float64, v::Float64, p::Float64, ε::Float64 )
        return new(b, c, β, u, v, p, ε, [0.0 b; -c b-c])
    end
end

mutable struct Simulation
    # structure storing information about the simulation
    # this is temporary storage *within* a generation
    # any info to be carried over for tracking should be in pop

    i_action::Int64                   # temporarily store action of i in a given pairwise game
    j_action::Int64                   # temporarily store action of j in a given pairwise game
    learner::Int64                    # temporarily store learner in a given generation
    non_learners::Array{Int64,1}      # temporarily store list of individuals other than the learner
    role_model::Int64                 # temporarily store role model in a given generation
    strategy_copy::Array{Int64,1}     # temporarily store learner's strategies before imitation
    sets_copy::Array{Int64,1}         # temporarily store learner's opinions before imitation
    fitnesses::Array{Float64,1}       # temporarily store fitnesses of everyone but the learner
    cityblock_ind::Array{Int64,2}     # temporarily store pairwise cityblock distance between individuals
    cityblock_party::Array{Int64,2}   # temporarily store within-party pairwise cityblock distance between individuals
    hamming_ind::Array{Int64,2}       # temporarily store pairwise hamming distance between individuals
    hamming_party::Array{Int64,2}     # temporarily store within-party hamming cityblock distance between individuals
    cooperation::Array{Int64, 3}      # temporarily store within-party or between-party cooperation

    # constructor
    function Simulation( sets::Sets )
        # initialize
        i_action, j_action, learner, role_model = -1, -1, -1, -1  # initializing with a value other than 0(D)/1(C)
        non_learners    = zeros(Int64, sets.N-1)
        strategy_copy   = [-1, -1]
        sets_copy       = -ones(Int64, sets.M)
        fitnesses       = -ones(Float64, sets.N-1)

        cityblock_ind   = zeros(Int64, sets.N, sets.N) 
        cityblock_party = zeros(Int64, sets.N÷2, sets.N÷2)
        hamming_ind     = zeros(Int64, sets.N, sets.N) 
        hamming_party   = zeros(Int64, sets.N÷2, sets.N÷2)
        cooperation     = zeros(Int64, sets.M, sets.N÷2, sets.N÷2)

        return new(i_action, j_action, learner, non_learners, role_model, strategy_copy, sets_copy, fitnesses,
                   cityblock_ind, cityblock_party, hamming_ind, hamming_party, cooperation)
    end
end

mutable struct Population
    # structure storing information about the population
    
    sets::Sets                        # the set distribution
    game::Game                        # not mutable: game attributes and parameters
    sim::Simulation                   # temporary info in simulations
    strategies::Array{Int64, 2}       # array of game strategies - 0:DD, 1:DC, 2:CD, 3:CC
    payoffs::Array{Float64, 1}        # array of payoffs
    all_strategies::Array{Int64, 1}   # which strategies are allowed to appear
    generation::Int64                 # current generation
    verbose::Bool                     # turn this on for error tracking
    sets_changed::Bool                # track changes in sets
    strategies_changed::Bool          # track changes in strategies
    prev_actions::Array{Int64, 3}     # the last action each individual took toward each other
    prev_learner::Int64               # track learner in the previous round
    payoffs_mat::Array{Float64, 2}    # track pairwise payoffs
    all_opinions::Array{Int64, 1}     # list of all possible opinions (in base10)
    all_memberships::Array{Int64,1}   # list of all possible set memberships (in base10)
    strategies_base10::Array{Int64,1} # list of strategies in base 10
    prev_interactions::Int64          # track number of interactions
    num_imitating::Int64              # track number of imitation events
    num_not_imitating::Int64          # track number of !(imitation events)
    num_strategies_changed::Int64     # track number of strategy changes
    num_sets_changed::Int64           # track number of opinion changes
    num_both_changed::Int64           # track number of times when both strategy + opinions change
    num_random_opinions::Int64        # track number of times random opinions were adopted
    num_biased_opinions::Int64        # track number of times biased opinions were adopted
    num_random_strategies::Int64      # track number of times random strategies were adopted

    # constructor for the case where initial strategy distribution is specified
    function Population( sets::Sets, game::Game, initial_strategies::Array{Int64, 2}, verbose::Bool=false, sets_changed::Bool=true,    strategies_changed::Bool=true )

        # initializing the population with specified strategies
        payoffs      = zeros(Float64, sets.N)
        generation   = 0
        payoffs_mat  = zeros(Float64, (sets.N, sets.N))
        
        # conditional: 00 (DD), 01 (DC), 01 (CD), or 11 (CC) (in binary)
        strategies        = initial_strategies
        all_strategies    = collect(0:3)

        # misc
        prev_actions      = zeros(Int64, sets.M, sets.N, sets.N)
        all_opinions      = (-3^sets.M+1)÷2 : (3^sets.M-1)÷2
        all_memberships   = 1 : (3^sets.M-1)÷2
        strategies_base10 = [2*strategies[i,1] + strategies[i,2] for i in 1:sets.N]
        prev_interactions, prev_learner = 0, 0

        # overall trackers
        num_imitating, num_not_imitating, num_strategies_changed, num_sets_changed = 0, 0, 0, 0
        num_random_opinions, num_biased_opinions, num_random_strategies            = 0, 0, 0

        # initialize sim::Simulation
        sim = Simulation(sets)

        return new(sets, game, sim, strategies, payoffs, all_strategies, 
            generation, verbose, sets_changed, strategies_changed, 
            prev_actions, prev_learner, payoffs_mat, 
            all_opinions, all_memberships, strategies_base10, prev_interactions,
            num_imitating, num_not_imitating, num_strategies_changed, num_sets_changed, 
            num_both_changed, num_random_opinions, num_biased_opinions, num_random_strategies)
    end

    # same constructor as above, but allows randomized initial strategies
    function Population( sets::Sets, game::Game, verbose::Bool=false, sets_changed::Bool=true, strategies_changed::Bool=true )

        # begin by initializing the population with random strategies
        payoffs      = zeros(Float64, sets.N)
        generation   = 0
        payoffs_mat  = zeros(Float64, (sets.N, sets.N))
        
        # conditional: 00 (DD), 01 (DC), 01 (CD), or 11 (CC) (in binary)
        strategies        = rand(collect(0:1), sets.N, 2)  # each strat is 0 or 1, total 4
        all_strategies    = collect(0:3)

        # misc
        prev_actions      = zeros(Int64, sets.M, sets.N, sets.N)
        all_opinions      = (-3^sets.M+1)÷2 : (3^sets.M-1)÷2
        all_memberships   = 1 : (3^sets.M-1)÷2
        strategies_base10 = [2*strategies[i,1] + strategies[i,2] for i in 1:sets.N]
        prev_interactions, prev_learner = 0, 0

        # overall trackers
        num_imitating, num_not_imitating, num_strategies_changed, num_sets_changed = 0, 0, 0, 0
        num_random_opinions, num_biased_opinions, num_random_strategies            = 0, 0, 0

        # initialize sim::Simulation
        sim = Simulation(sets)

        return new(sets, game, sim, strategies, payoffs, all_strategies, 
            generation, verbose, sets_changed, strategies_changed, 
            prev_actions, prev_learner, payoffs_mat, 
            all_opinions, all_memberships, strategies_base10, prev_interactions,
            num_imitating, num_not_imitating, num_strategies_changed, num_sets_changed, 
            num_both_changed, num_random_opinions, num_biased_opinions, num_random_strategies)
    end
end
