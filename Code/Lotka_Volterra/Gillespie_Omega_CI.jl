
## Load Packages
using Omega
using StatsBase
using Random
using Plots
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

plot_vals = hcat(values(samples[2][6]["prey"]), values(samples[2][6]["pred"]))
plot(1:10001, plot_vals[1:10001,1:2],
        title = "Basic Simulation",
        label = ["Prey" "Pred"],
        lw = 1.25)

## Sample with intervention
samples_int = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, gillespie),
                prey_init > 120, 5, alg = RejectionSample)

plot_vals = hcat(values(samples[2][6]["prey"]), values(samples_int[2][6]["prey"]))
plot(1:10001, plot_vals[1:10001,1:2],
        title = "Intervention Comparison",
        label = ["Prey No Intervention" "Prey > 120"],
        lw = 1.25)

plot_vals = hcat(values(samples_int[2][6]["prey"]), values(samples_int[2][6]["pred"]))
plot(1:10001, plot_vals[1:10001,1:2],
        title = "Intervention Simulation",
        label = ["Prey" "Pred"],
        lw = 1.25)

# Histogram
function intervention_hist(sim1, sim2, b)

    sim1_final_vals = zeros(0)
    sim2_final_vals = zeros(0)
    for i in 1:length(sim1)
        append!(sim1_final_vals, sim1[1][6]["prey"][i])
        append!(sim2_final_vals, sim2[1][6]["prey"][i])
    end

    final_val_diff = sim1_final_vals - sim2_final_vals
    histogram(final_val_diff, bins = b)
end

samples_100 = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, gillespie),
                100, alg = RejectionSample)
samples_int_100 = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, gillespie),
                prey_init > 120, 100, alg = RejectionSample)

intervention_hist(samples,samples_int, 20)

# Sample Rate
samples_rate_int = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, gillespie),
                    spawn_prey > 1.2, 5, alg = RejectionSample)

plot_vals = hcat(values(samples[2][6]["prey"]), values(samples_rate_int[2][6]["prey"]))
plot(1:5001, plot_vals[1:5001,1:2],
        title = "Intervention Comparison",
        label = ["Prey No Intervention" "Prey - spawn_prey_rate > 1.2"],
        lw = 1.25)

## Idea for intervention at t_now
function sample_weights(rng, hazards)
    result = sample(collect(keys(hazards)), Weights(collect(values(hazards))))
    return result
end

function initial_simulation(rng, t, transitions)

    ecology = Dict("prey" => prey_init(rng), "pred" => pred_init(rng))
    theta = Dict("spawn_prey" => spawn_prey(rng),
                 "prey2pred" => prey2pred(rng),
                 "pred_dies" => pred_dies(rng))

    hazards = get_hazards(ecology, theta)

    # Does this make sense?
    hazard_sample = ciid(sample_weights, hazards)
    transition = transitions[hazard_sample(rng)]
    t_new = t + sum(values(hazards))

    ecology_new_prey = ecology["prey"] + transition[1]
    ecology_new_pred = ecology["pred"] + transition[2]
    # Enforce only positive integers
    ecology_new_prey = max(1, ecology_new_prey)
    ecology_new_pred = max(1, ecology_new_pred)

    return Dict("prey" => ecology_new_prey,
                "pred" => ecology_new_pred,
                "t" => t_new)
end

temp = ciid(initial_simulation, t0, transitions)
samples = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, temp),
                5, alg = RejectionSample)

function get_hazards_fix(ecology, theta)
    """
    A function to fix annoying julia things
    """
    return Dict(
        "spawn_prey" => theta["spawn_prey"][1] * ecology["prey"],
        "prey2pred" => theta["prey2pred"][1] * ecology["prey"] * ecology["pred"],
        "pred_dies" => theta["pred_dies"][1] * ecology["pred"]
        )
end

function one_simulation(rng, n, theta, transitions)
    ecology = functions_list[n](rng)
    ecology_temp = copy(ecology)
    t = ecology["t"]
    delete!(ecology_temp, "t")
    hazards = get_hazards_fix(ecology_temp, theta)
    # Does this make sense?
    #hazard_sample = ciid(sample_weights, hazards)
    transition = transitions[sample(collect(keys(hazards)), Weights(collect(values(hazards))))]
    t_new = t + sum(values(hazards))

    ecology_new_prey = ecology_temp["prey"] + transition[1]
    ecology_new_pred = ecology_temp["pred"] + transition[2]
    # Enforce only positive integers
    ecology_new_prey = max(1, ecology_new_prey)
    ecology_new_pred = max(1, ecology_new_pred)

    return Dict("prey" => ecology_new_prey,
                "pred" => ecology_new_pred,
                "t" => t_new)
end

temp = ciid(initial_simulation, t0, transitions)
functions_list = Any[]
push!(functions_list, temp)
N = 5
for f in 2:N
    last = f - 1
    temp = ciid(one_simulation, last, theta, transitions)
    push!(functions_list, temp)
end

push!(functions_list,prey_init)
push!(functions_list, pred_init)
push!(functions_list, spawn_prey)
push!(functions_list, prey2pred)
push!(functions_list, pred_dies)

samples = rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies,
functions_list[1], functions_list[2]),
                5, alg = RejectionSample)

samples = rand(Tuple(x for x in functions_list),
                5, alg = RejectionSample)
