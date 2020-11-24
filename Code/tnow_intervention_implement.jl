## Load Packages
using Omega
using StatsBase
using Random
using Plots
using Distributions

Random.seed!(123)

## Define functions
function sample_weights(rng, hazards)
    result = sample(collect(keys(hazards)), Weights(collect(values(hazards))))
    return result
end

function initial_simulation(rng, t, transitions)

    """
    Simulates first run of Gillespie. Needs to be different because
    it takes in the intial prey, pred and rates.
    """

    ecology = Dict("prey" => prey_init(rng), "pred" => pred_init(rng))
    # theta = Dict("spawn_prey" => spawn_prey(rng),
    #              "prey2pred" => prey2pred(rng),
    #              "pred_dies" => pred_dies(rng))

    hazards = get_hazards_fix(ecology, theta)

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

function one_simulation(rng, n, theta, transitions)

    """
    Simulates all steps of gillespie except intial step. Takes into account
    the results of the previous step "n".
    """

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

function get_hazards_fix(ecology, theta)
    """
    Same as above, for some reason rates become lists in subsequent runs so
    this fixes that. Juliaaaa issues 😫
    """
    return Dict(
        "spawn_prey" => theta["spawn_prey"][1] * ecology["prey"],
        "prey2pred" => theta["prey2pred"][1] * ecology["prey"] * ecology["pred"],
        "pred_dies" => theta["pred_dies"][1] * ecology["pred"]
        )
end

## Run simulation
# Define gillespie parameters
Pre = [[1, 0], [1, 1], [0, 1]]
Post = [[2, 0], [0, 2], [0, 0]]
transition_mat = Post - Pre
transitions = Dict("spawn_prey" => transition_mat[1,],
                    "prey2pred" => transition_mat[2,],
                    "pred_dies" => transition_mat[3,])
t0=0.0
N = 10000

# Random variables for starting values
prey_init = normal(100, 10)
pred_init = normal(50, 10)

# Random variables for rates
spawn_prey = normal(1, .1)
prey2pred = normal(.01, .001)
pred_dies = normal(.5, .05)

theta = Dict("spawn_prey" => rand(spawn_prey,1),
             "prey2pred" => rand(prey2pred,1),
             "pred_dies" => rand(pred_dies,1))

# Random variables for each step of simulation
temp = ciid(initial_simulation, t0, transitions) # initial step
functions_list = Any[]
push!(functions_list, temp)

for f in 2:N
    last = f - 1
    temp = ciid(one_simulation, last, theta, transitions) # individual step
    push!(functions_list, temp)
end

# Add initialized random vars to list
push!(functions_list,prey_init)
push!(functions_list, pred_init)
push!(functions_list, spawn_prey)
push!(functions_list, prey2pred)
push!(functions_list, pred_dies)

# sample
samples = rand(Tuple(x for x in functions_list),
                5, alg = RejectionSample)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:N
    push!(prey_vals,samples[3][x]["prey"])
    push!(pred_vals,samples[3][x]["pred"])
end

plot(hcat(prey_vals,pred_vals))

# sample
samples = rand(Tuple(x for x in functions_list),prey_init > 120
                ,, alg = RejectionSample)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:N
    push!(prey_vals,samples[3][x]["prey"])
    push!(pred_vals,samples[3][x]["pred"])
end

plot(hcat(prey_vals,pred_vals))
