# plot.jl: functions for plotting
using PyCall: PyObject, PyDict
using PyPlot: matplotlib, Figure, plt, withfig

# plotting parameters
rcParams = PyDict(matplotlib."rcParams")
rcParams["font.sans-serif"] = ["Arial"];
rcParams["font.size"] = 12;
colpalette = ["red", "orange", "lightblue", "blue"];

# functions for plotting
function plot_strat_dist(
    pop::Population,                        # need pop.sets.conditional
    data::Array{Float64, 2},                # data to plot
    ax::PyObject,                           # axis on which to plot
    ylabel::String="Frequency",             # y-axis label # change for cumulative distribution
    gen_start::Int64=1,                     # starting time step
    gen_skip::Int64=1,                      # time step interval
    generations::Int64=size(data,1),        # final time step / total sim length
    xlabel::String="Time step",             # x-axis label
    colpalette::Array{String, 1}=colpalette # color palette
    )
    # function to plot strategy distribution 
    strat_ids = ["DD","DC","CD","CC"] # names of strategies

    [ax.plot(gen_start:gen_skip:generations, data[gen_start:gen_skip:generations,x], 
        label=strat_ids[x], linewidth=0.75, 
        color=colpalette[x]) for x in 1:length(pop.all_strategies)]
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim([gen_start,generations])
    ax.set_ylim([-0.05,1.05])
    ax.legend(bbox_to_anchor=(1.15, 0.5), loc=5, borderaxespad=0., fontsize=10)

end

# another method, with effective coopeartion
function plot_strat_dist(
    pop::Population,                        # need pop.sets.conditional
    data::Array{Float64, 2},                # data to plot
    data2::Array{Float64, 1},               # data to plot
    ax::PyObject,                           # axis on which to plot
    ylabel::String="Frequency",             # y-axis label # change for cumulative distribution
    gen_start::Int64=1,                     # starting time step
    gen_skip::Int64=1,                      # time step interval
    generations::Int64=size(data,1),        # final time step / total sim length
    xlabel::String="Time step",             # x-axis label
    colpalette::Array{String, 1}=colpalette # color palette
    )
    # function to plot strategy distribution 
    strat_ids = ["DD","DC","CD","CC","Cooperation"] # names of strategies

    [ax.plot(gen_start:gen_skip:generations, data[gen_start:gen_skip:generations,x], 
        label=strat_ids[x], linewidth=0.75, 
        color=colpalette[x]) for x in 1:length(pop.all_strategies)]
    ax.plot(gen_start:gen_skip:generations, data2[gen_start:gen_skip:generations], 
        label=strat_ids[5], linewidth=0.5, color="Black", linestyle="dotted")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim([gen_start,generations])
    ax.set_ylim([-0.05,1.05])
    ax.legend(bbox_to_anchor=(1.2, 0.5), loc=5, borderaxespad=0., fontsize=8)

end

function plot_set_dist(
    pop::Population,                     # need pop.sets.conditional
    data::Array{Int64, 2},               # data to plot
    ax::PyObject,                        # axis on which to plot
    ylabel::String="No. of individuals", # y-axis label # change for cumulative distribution
    gen_start::Int64=1,                  # starting time step
    gen_skip::Int64=1,                   # time step interval
    generations::Int64=size(data,1),     # final time step / total sim length
    xlabel::String="Time step"           # x-axis label
    )
    # function to plot strategy distribution
    set_ids = ["Set $(x)" for x in 1:pop.sets.M] # define names of sets for legend

    [ax.plot(gen_start:gen_skip:generations, data[gen_start:gen_skip:generations,x], 
        label=set_ids[x], linewidth=0.75) for x in 1:(pop.sets.M)]
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim([gen_start,generations])
    ax.set_ylim([-1,(pop.sets.N)+1])
    ax.legend(bbox_to_anchor=(1.15, 0.5), loc=5, borderaxespad=0., fontsize=8) 

end

function plot_opn_dist(
    pop::Population,                          # need pop.sets.conditional
    data::Array{Float64, 2},                  # data to plot
    ax::PyObject,                             # axis on which to plot
    ylabel::String="Frequency",               # y-axis label # change for cumulative distribution
    gen_start::Int64=1,                       # starting time step
    gen_skip::Int64=1,                        # time step interval
    generations::Int64=size(data,1),          # final time step / total sim length
    xlabel::String="Time step"                # x-axis label
    )
    # function to plot average and variance of total opinions 
    opinion_ids::Array{String, 1} = ["$(base3todec(x, pop))" for x in pop.all_opinions]

    # ax.plot(gen_start:gen_skip:generations, data[gen_start:gen_skip:generations], linewidth=0.25)
    [ax.plot(gen_start:gen_skip:generations, 
        data[gen_start:gen_skip:generations,x], 
        label=opinion_ids[x], linewidth=0.75) for x in 1:length(pop.all_opinions)]
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim([gen_start,generations])
    ax.set_ylim([-0.05,1.05])
    ax.legend(bbox_to_anchor=(1.15, 0.5), loc=5, borderaxespad=0., fontsize=7) # note fontsize
end

function plot_total_opn(
    pop::Population,                          # need pop.sets.conditional
    data::Array{Float64, 1},                  # data to plot
    ax::PyObject,                             # axis on which to plot
    ylabel::String="Average\ntotal opinion",  # y-axis label # change for cumulative distribution
    gen_start::Int64=1,                       # starting time step
    gen_skip::Int64=1,                        # time step interval
    generations::Int64=size(data,1),          # final time step / total sim length
    xlabel::String="Time step"                # x-axis label
    )
    # function to plot average and variance of total opinions 
    ax.plot(gen_start:gen_skip:generations, 
            data[gen_start:gen_skip:generations], 
            label="value", linewidth=0.75, color="purple")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim([gen_start,generations])
    ax.set_ylim([min(0,minimum(data))-0.1, max(pop.sets.K,maximum(data))+0.1])
    ax.legend(bbox_to_anchor=(1.15, 0.5), loc=5, borderaxespad=0., fontsize=8)
end

function base3todec( base3::Int64, pop::Population )
    # convert opinion in int to vector of 1, 0, -1
    reverse( digits( base3 + vectobase3( repeat([1],pop.sets.M) ), base = 3, pad = pop.sets.M )) .- 1
end

function plotter( tracker::Tracker, gen_start::Int64, gen_skip::Int64, gen_end::Int64 )
    # function to make a quick plot
    fig, axs = plt.subplots(6, 2, figsize = (16,12)); 

    plot_strat_dist( tracker.pop, tracker.ens_strat_dist, tracker.ens_total_cooperation, axs[1], "Frequency", gen_start, gen_skip, gen_end ); 
    plot_opn_dist( tracker.pop, tracker.ens_opn_dist, axs[3], "Frequency", gen_start, gen_skip, gen_end); # opinion distribution
    plot_set_dist( tracker.pop, tracker.ens_pop_dist, axs[9], "No. of individuals", gen_start, gen_skip, gen_end );
    plot_strat_dist( tracker.pop, tracker.ens_running_strat_dist, axs[2], "Running\navg frequency", gen_start, gen_skip, gen_end ); 
    plot_total_opn( tracker.pop, tracker.ens_total_means, axs[7], "Average\ntotal opinion", gen_start, gen_skip, gen_end );  
    plot_total_opn( tracker.pop, tracker.ens_total_vars, axs[8], "Variance\ntotal opinion", gen_start, gen_skip, gen_end );
    # plot_total_opn( tracker.pop, tracker.ens_total_cooperation, axs[7], "Effective cooperation", gen_start, gen_skip, gen_end );
    plot_total_opn( tracker.pop, tracker.ens_opn_simpson, axs[4], "Opinion diversity\n(1 = low)", gen_start, gen_skip, gen_end );
    plot_total_opn( tracker.pop, tracker.ens_sets_simpson, axs[10], "Set diversity\n (1 = low)", gen_start, gen_skip, gen_end );
    plot_total_opn( tracker.pop, tracker.ens_hamming_ind, axs[5], "Pairwise\nHamming distance", gen_start, gen_skip, gen_end );
    plot_total_opn( tracker.pop, tracker.ens_hamming_party, axs[6], "Within-party\nHamming distance", 
    gen_start, gen_skip, gen_end );
    plot_total_opn( tracker.pop, tracker.ens_cityblock_ind, axs[11], "Pairwise\nCityblock distance", gen_start, gen_skip, gen_end );
    plot_total_opn( tracker.pop, tracker.ens_cityblock_party, axs[12], "Within-party\nCityblock distance", gen_start, gen_skip, gen_end );

    return fig 
end