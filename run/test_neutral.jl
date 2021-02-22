using CooperationPolarization
using BenchmarkTools
using Profile
using PyPlot: plt, savefig
using Dates

const b = 1.0 # benefit to cooperating
const c = 0.2 # cost to cooperating

const N = 40 # population size
const M = 1 # number of sets
const K = 1 # number of issues to care about
const P = 2 # number of parties

const β = 0.000    # selection strength
const u = 0.001    # strategy mutation rate
const v = 0.001    # set mutation rate
const p = 0.0     # party bias (0 = low, 1 = high)
const ϵ = 1.       # damping factor (0 = no bias in mutation, 1 = max bias (= p) in mutation)
# note: α = 0 automatically

# simulation parameters
const verbose       = false
const verbose_track = false
const generations   = 10000000
const testtype      = "test_neutral"

# plotting parameters
const gen_skip  = 1 # max( div(generations, 10000), 1)
const gen_start = max( 1, generations-1000000 )
const gen_end   = generations

const start_time = Dates.format(now(),"yymmdd-HHMMSS.sss")

println("parameters defined")

sets      = random_sets(N, M, K, P)
game      = Game(b, c, β, u, v, p, ϵ)
pop       = Population(sets, game, verbose)
n_tracker = NeutralTracker(pop, generations, verbose_track)

println("Payoff matrix: \t $(game.A)")
println("Payoffs:\t $(pop.payoffs)")
println("Verbose?:\t $(pop.verbose)")
println("Initial:\t strategies_changed = $(pop.strategies_changed), sets_changed = $(pop.sets_changed)")

@time for gen::Int64 in 1:generations
    if generations > 10 && (gen-1) % (generations ÷ 10) == 0  # set up to print 100 generations
        println("\tinitiating generation $(n_tracker.pop.generation)")
    end
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

# plotting
# fig = plotter( n_tracker, gen_start, gen_skip, gen_end )

fig, axs = plt.subplots(4, 4, figsize = (15,8))

plot_total_opn( n_tracker.pop, n_tracker.ens_y_running, axs[1], "y, runnig avg", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_z_running, axs[2], "z, runnig avg", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_g_running, axs[3], "g, runnig avg", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_h_running, axs[4], "h, runnig avg", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_sia_sjd_running, axs[9],  "<s_ia s_jd>, runnig avg", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_sia_sid_running, axs[10], "<s_ia s_id>, runnig avg", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_y, axs[5], "y", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_z, axs[6], "z", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_g, axs[7], "g", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_h, axs[8], "h", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_sia_sjd, axs[13], "<s_ia s_jd>", gen_start, gen_skip, gen_end)
plot_total_opn( n_tracker.pop, n_tracker.ens_sia_sid, axs[14], "<s_ia s_id>", gen_start, gen_skip, gen_end)

plot_strat_dist( n_tracker.pop, n_tracker.ens_strat_dist, axs[15], "abundance", gen_start, gen_skip, gen_end)
plot_strat_dist( n_tracker.pop, n_tracker.ens_strat_running, axs[11], "abundance", gen_start, gen_skip, gen_end)

plt.suptitle("Parameters: b=$b, c=$c, N=$N, M=$M, K=$K, β=$β, u=$u, v=$v, p=$p, ϵ=$ϵ", fontsize=12)
fig.tight_layout(rect=[0, 0.03, 1, 0.97])

display(fig)

# save file
# savefig("plots/$(testtype)_$(filename).pdf")
