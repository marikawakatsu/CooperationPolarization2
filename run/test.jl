using CooperationPolarization
using BenchmarkTools
using Profile
using PyPlot: plt, savefig
using Dates

const b = 200. # benefit to cooperating
const c = 0.2 # cost to cooperating

const N = 40 # population size
const M = 1 # number of sets
const K = 1 # number of issues to care about
const P = 2 # number of parties

const β = 0.001    # selection strength
const u = 0.001    # strategy mutation rate
const v = 0.001    # set mutation rate
const p = 1.       # party bias (0 = low, 1 = high)
const ϵ = 1.       # damping factor (0 = no bias in mutation, 1 = max bias (= p) in mutation)
# note: α = 0 automatically

# simulation parameters
const verbose       = false
const verbose_track = false
const generations   = 10000000
const testtype      = "test"

# plotting parameters
const gen_skip  = 1 # max( div(generations, 10000), 1)
const gen_start = max( 1, generations-1000000 )
const gen_end   = generations

const start_time = Dates.format(now(),"yymmdd-HHMMSS.sss")

println("parameters defined")

sets    = random_sets(N, M, K, P)
game    = Game(b, c, β, u, v, p, ϵ)
pop     = Population(sets, game, verbose)
tracker = Tracker(pop, generations, verbose_track)

println("Payoff matrix: \t $(game.A)")
println("Payoffs:\t $(pop.payoffs)")
println("Verbose?:\t $(pop.verbose)")
println("Initial:\t strategies_changed = $(pop.strategies_changed), sets_changed = $(pop.sets_changed)")

@time for g::Int64 in 1:generations
    if generations > 10 && (g-1) % (generations ÷ 10) == 0  # set up to print 100 generations
        println("\tinitiating generation $(tracker.pop.generation)")
    end
    # update population
    evolve!(tracker.pop, 1)  # this is the same as evolve!(pop, 1)
    track!(tracker, g)

    # check pop and tracker.pop are the same
    # println("\t$(pop.sets.h == tracker.pop.sets.h)") # just to check

    # check whether the updating step is working properly
    # if true, there should be a mix of 'true' and 'false' for sets_chagned and strategies_changed
    # println("\tsets_changed = $(pop.sets_changed), \tstrat_changed = $(pop.strategies_changed)") 

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

# plotting
fig = plotter( tracker, gen_start, gen_skip, gen_end )

plt.suptitle("Parameters: b=$b, c=$c, N=$N, M=$M, K=$K, β=$β, u=$u, v=$v, p=$p, ϵ=$ϵ", fontsize=12)
fig.tight_layout(rect=[0, 0.03, 1, 0.97])

display(fig)

# save file
# savefig("plots/$(testtype)_$(filename).pdf")
