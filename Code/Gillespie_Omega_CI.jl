
## Load Packages
using Omega
using UnicodePlots
using StatsBase
using Plots

## Define Functions
function gillespie(rng)

    """
    Implementation of the Gillespie algorithm for Lotka-Volterra model
    """
    Pre = [[1, 0], [1, 1], [0, 1]]
    Post = [[2, 0], [0, 2], [0, 0]]
    transition_mat = Post - Pre
    transitions = Dict("spawn_prey" => transition_mat[1,],
                        "prey2pred" => transition_mat[2,],
                        "pred_dies" => transition_mat[3,])
    t0=0.0
    N = 10000
    initial = Dict("prey" => prey_init(rng), "pred" => pred_init(rng))
    trajectory = Dict("prey" => [initial["prey"]], "pred" => [initial["pred"]])
    times = [t0]
    ecology = initial

    theta = Dict("spawn_prey" => spawn_prey(rng),
                 "prey2pred" => prey2pred(rng),
                 "pred_dies" => pred_dies(rng))
    t = times[1]

    for i = 1:N
        hazards = get_hazards(ecology, theta)
        transition = transitions[sample(collect(keys(hazards)), Weights(collect(values(hazards))))]
        t = t + sum(values(hazards))

        ecology["prey"] += transition[1]
        ecology["pred"] += transition[2]
        # Enforce only positive integers
        ecology["prey"] = max(0, ecology["prey"])
        ecology["pred"] = max(0, ecology["pred"])

        append!(trajectory["prey"], ecology["prey"])
        append!(trajectory["pred"], ecology["pred"])
        append!(times, t)
    end

    return Dict(
        "prey" => trajectory["prey"],
        "pred" => trajectory["pred"],
        "times" => times
        )
end

function get_hazards(ecology, theta)
    """
    Compute the hazard function given the current states.  Note that in this
        model the function depends only on state, not time, though this is not
        the case in general. "spawn_prey" represents the event of a prey being born,
        "prey2pred" represents a predator consuming a new prey and consequently spawning
        a new predator, "pred_dies" represents the death of a predator.

    args:
        ecology(dict): a dictionary containing species counts.
        theta(dict): A dictionary where keys represent events (reactions) and
            values are the location hyperparameters for rate parameter constants
            corresponding to those events.
    """

    return Dict(
        "spawn_prey" => theta["spawn_prey"] * ecology["prey"],
        "prey2pred" => theta["prey2pred"] * ecology["prey"] * ecology["pred"],
        "pred_dies" => theta["pred_dies"] * ecology["pred"]
        )
end

# Initiate starting values randomely
prey_init = normal(100, 10)
pred_init = normal(50, 10)

# Should these be random variables? maybe not? they are in the paper
spawn_prey = uniform([1, .9, 1.5])
prey2pred = uniform([.005, .004, .01])
pred_dies = uniform([.6, .4, .8])

gillespie_ = ciid(gillespie)
samples = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, gillespie_),
                5, alg = RejectionSample)

samples_int = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, gillespie_),
                    prey_init > 125, 5, alg = RejectionSample)

plot_vals = hcat(values(samples[2][6]["prey"]), values(samples_int[2][6]["prey"]))
plot(1:10001, plot_vals)

plot_vals = hcat(values(samples[1][6]["prey"]), values(samples[1][6]["pred"]))
plot(1:10001, plot_vals)
