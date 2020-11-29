## Load Packages
using Omega
using StatsBase
using Random
using Plots
using Distributions

Random.seed!(123)

function one_simulation_prey(rng, n, transitions)

    """
    Simulates one step of gillespie for prey.
    """

    #t = t_list[n](rng)
    hazard_result = hazards_list[n](rng)
    prey_val = prey_list[n](rng)
    transition = transitions[hazard_result]
    new_prey = prey_val + transition[1]

    # Enforce only positive integers
    new_prey = max(1, new_prey)
    return new_prey
end

function one_simulation_pred(rng, n, transitions)

    """
    Simulates one step of gillespie for pred.
    """

    #t = t_list[n](rng)
    hazard_result = hazards_list[n](rng)
    pred_val = pred_list[n](rng)
    transition = transitions[hazard_result]
    new_pred = pred_val + transition[2]

    # Enforce only positive integers
    new_pred = max(1, new_pred)
    return new_pred
end

function get_hazards(rng, n, theta)
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

    ecology = Dict("prey" => prey_list[n](rng), "pred" => pred_list[n](rng))

    hazards = Dict(
        "spawn_prey" => theta["spawn_prey"][1] * ecology["prey"],
        "prey2pred" => theta["prey2pred"][1] * ecology["prey"] * ecology["pred"],
        "pred_dies" => theta["pred_dies"][1] * ecology["pred"]
        )
    result = sample(collect(keys(hazards)), Weights(collect(values(hazards))))
    return result
end


## Run simulation
# Define gillespie parameters
Pre = [[1, 0], [1, 1], [0, 1]]
Post = [[2, 0], [0, 2], [0, 0]]
transition_mat = Post - Pre
transitions = Dict("spawn_prey" => transition_mat[1,],
                    "prey2pred" => transition_mat[2,],
                    "pred_dies" => transition_mat[3,])
#t0=0.0
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

# Function to join lists together
@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = tuplejoin(tuplejoin(x, y), z...)

## Sample without intervention
# Random variables for each step of simulation
hazards_list = Any[]
prey_list = Any[]
pred_list = Any[]
push!(prey_list, prey_init)
push!(pred_list, pred_init)

for f in 2:N
    last = f - 1
    hazards_temp = ciid(get_hazards, last, theta) # individual step
    prey_temp = ciid(one_simulation_prey, last, transitions) # individual step
    pred_temp = ciid(one_simulation_pred, last, transitions) # individual step
    push!(hazards_list, hazards_temp)
    push!(prey_list, prey_temp)
    push!(pred_list, pred_temp)
end

samples = rand(tuplejoin(Tuple(x for x in hazards_list),
                Tuple(x for x in prey_list),
                Tuple(x for x in pred_list)),
                5, alg = RejectionSample)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,samples[2][N+x])
    push!(pred_vals,samples[2][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Intervention Simulation",
        label = ["Prey" "Pred"],
        lw = 1.25)

## Sample with intervention
hazards_list = Any[]
prey_list = Any[]
pred_list = Any[]
push!(prey_list, prey_init)
push!(pred_list, pred_init)

for f in 2:N
    last = f - 1
    hazards_temp = ciid(get_hazards, last, theta) # individual step
    prey_temp = ciid(one_simulation_prey, last, transitions) # individual step
    pred_temp = ciid(one_simulation_pred, last, transitions) # individual step
    push!(hazards_list, hazards_temp)
    push!(prey_list, prey_temp)
    push!(pred_list, pred_temp)
end

random_var_tuple = tuplejoin(Tuple(x for x in hazards_list),
                Tuple(x for x in prey_list),
                Tuple(x for x in pred_list))
samples = rand(random_var_tuple, random_var_tuple[15000] > 400.0,
                100, alg = RejectionSample)

compile = []
for x in 1:100
    push!(compile, samples[x][15000])
end
histogram(compile)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,samples[3][N+x])
    push!(pred_vals,samples[3][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Intervention Simulation",
        label = ["Prey" "Pred"],
        lw = 1.25)
