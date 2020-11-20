
## Load Packages
using Omega
using UnicodePlots
using StatsBase
using Plots
using Random
using Distributions

Random.seed!(123)

## Define Functions
function gillespie_(rng, transitions, N, t0=0.0)

    """
    Implementation of the Gillespie algorithm for Lotka-Volterra model

    Args:
    transitions(dict): A dictionary with keys "prey" and "pred", and each element
        is a list containing the row of the transition matrix.
    N(int): Number of desired simulation iterations.
    t0(float): initial time point.
    """

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
        # Can we replace this with a random variable?
        # It already samples one of three hazards according to some weights
        transition = transitions[sample(collect(keys(hazards)), Weights(collect(values(hazards))))]
        t = t + sum(values(hazards))

        ecology["prey"] += transition[1]
        ecology["pred"] += transition[2]
        # Enforce only positive integers
        ecology["prey"] = max(1, ecology["prey"])
        ecology["pred"] = max(1, ecology["pred"])

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

# Define gillespie parameters
Pre = [[1, 0], [1, 1], [0, 1]]
Post = [[2, 0], [0, 2], [0, 0]]
transition_mat = Post - Pre
transitions = Dict("spawn_prey" => transition_mat[1,],
                    "prey2pred" => transition_mat[2,],
                    "pred_dies" => transition_mat[3,])
t0=0.0
N = 10000

# Initiate starting values randomely
prey_init = normal(100, 10)
pred_init = normal(50, 10)

# Should these be random variables? maybe not? they are in the paper
spawn_prey = normal(1, .1)
prey2pred = normal(.01, .001)
pred_dies = normal(.5, .05)

# Make gillespie into random variable
gillespie = ciid(gillespie_, transitions, 10000, 0.0)

# Sample from random variables
samples = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, gillespie),
                5, alg = RejectionSample)


plot_vals = hcat(values(samples[1][6]["prey"]), values(samples[1][6]["pred"]))
plot(1:5001, plot_vals[1:5001,1:2],
        title = "Basic Simulation",
        label = ["Prey" "Pred"],
        lw = 1.25)

## Sample with intervention
samples_int = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, gillespie),
                    prey_init > 125, 5, alg = RejectionSample)

plot_vals = hcat(values(samples[2][6]["prey"]), values(samples_int[2][6]["prey"]))
plot(1:5001, plot_vals[1:5001,1:2],
        title = "Intervention Comparison",
        label = ["Prey No Intervention" "Prey > 125"],
        lw = 1.25)

samples_rate_int = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, gillespie),
                    spawn_prey > 1.2, 5, alg = RejectionSample)

plot_vals = hcat(values(samples[2][6]["prey"]), values(samples_rate_int[2][6]["prey"]))
plot(1:5001, plot_vals[1:5001,1:2],
        title = "Intervention Comparison",
        label = ["Prey No Intervention" "Prey - spawn_prey_rate > 1.2"],
        lw = 1.25)

## Counterfactual
gillespie_new = replace(gillespie, prey_init => 120)
counter_samples = rand(gillespie_new, 5, alg = RejectionSample)

plot_vals = hcat(values(samples[2][6]["prey"]), values(counter_samples[2]["prey"]))
plot(1:5001, plot_vals[1:5001,1:2],
        title = "Counterfactual Comparison",
        label = ["Prey No Intervention" "Prey > 125"],
        lw = 1.25)

## Idea for intervention at t_now
#   Create a new gillespie function (gillespie_new) that
#   has a paramter for changing a value at t_now. Use replace
#   to update simulation from old gillespie to gillespie_new.
#   This should give all the same values until t_now, which
#   then will update to whatever intervention was set.
#   I think this is how they do it in the Omega paper.
function gillespie_time_intervention_(rng, transitions, N, time, update,t0=0.0)

    """
    Implementation of the Gillespie algorithm for Lotka-Volterra model

    Args:
    transitions(dict): A dictionary with keys "prey" and "pred", and each element
        is a list containing the row of the transition matrix.
    N(int): Number of desired simulation iterations.
    t0(float): initial time point.
    """

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
        # Can we replace this with a random variable?
        # It already samples one of three hazards according to some weights
        transition = transitions[sample(collect(keys(hazards)), Weights(collect(values(hazards))))]
        t = t + sum(values(hazards))

        if N == time
            ecology["prey"] = update
            ecology["pred"] += transition[2]
        else
            ecology["prey"] += transition[1]
            ecology["pred"] += transition[2]
            # Enforce only positive integers
            ecology["prey"] = max(1, ecology["prey"])
            ecology["pred"] = max(1, ecology["pred"])
        end
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

gillespie_int = ciid(gillespie_time_intervention_, transitions, 10000, 2000, 100, 0.0)
gillespie_new = replace(gillespie, gillespie => gillespie_int)
counter_samples = rand(gillespie_new, 5, alg = RejectionSample)

plot_vals = hcat(values(samples[1][6]["prey"]), values(counter_samples[1]["prey"]))
plot(1:5001, plot_vals[1:5001,1:2],
        title = "Counterfactual Comparison",
        label = ["Prey No Intervention" "Prey > 125"],
        lw = 1.25)
