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
    labels = ["spawn_prey", "pred_dies","prey2pred"]
    transition = transitions[labels[hazard_result]]
    new_prey = prey_val + transition[1]

    # Enforce only positive integers
    max(1, new_prey)

end

function one_simulation_pred(rng, n, transitions)

    """
    Simulates one step of gillespie for pred.
    """

    #t = t_list[n](rng)
    hazard_result = hazards_list[n](rng)
    pred_val = pred_list[n](rng)
    labels = ["spawn_prey", "pred_dies","prey2pred"]
    transition = transitions[labels[hazard_result]]
    new_pred = pred_val + transition[2]

    # Enforce only positive integers
    max(1, new_pred)

end

function get_hazards(rng, n)
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
        "spawn_prey" => theta(rng)["spawn_prey"] * ecology["prey"],
        "prey2pred" => theta(rng)["prey2pred"] * ecology["prey"] * ecology["pred"],
        "pred_dies" => theta(rng)["pred_dies"] * ecology["pred"]
        )

    vals = collect(values(hazards))
    sum_vals = sum(vals)
    prob_vals = vals/sum_vals
    categorical(rng, prob_vals)
end

function generate_rates(rng)

    Dict("spawn_prey" => spawn_prey(rng),
         "prey2pred" => prey2pred(rng),
         "pred_dies" => pred_dies(rng))

end

## Initialize paramters
Pre = [[1, 0], [1, 1], [0, 1]]
Post = [[2, 0], [0, 2], [0, 0]]
transition_mat = Post - Pre
transitions = Dict("spawn_prey" => transition_mat[1,],
                    "prey2pred" => transition_mat[2,],
                    "pred_dies" => transition_mat[3,])
#t0=0.0
N = 10

# Random variables for starting values
prey_init = normal(100, 10)
pred_init = normal(50, 10)

# Random variables for rates
spawn_prey = normal(1., .1)
prey2pred = normal(.01, .001)
pred_dies = normal(.5, .05)

## Generate Model
hazards_list = Any[]
prey_list = Any[]
pred_list = Any[]
theta = ciid(generate_rates)
push!(prey_list, prey_init)
push!(pred_list, pred_init)
N = 10

for f in 2:N
    last = f - 1
    hazards_temp = ciid(get_hazards, last) # individual step
    prey_temp = ciid(one_simulation_prey, last, transitions) # individual step
    pred_temp = ciid(one_simulation_pred, last, transitions) # individual step
    push!(hazards_list, hazards_temp)
    push!(prey_list, prey_temp)
    push!(pred_list, pred_temp)
end

random_var_tuple = (Tuple(x for x in hazards_list)...,
                Tuple(x for x in prey_list)...,
                Tuple(x for x in pred_list)...,
                Tuple(Any[spawn_prey,prey2pred,pred_dies,theta])...)

## Sample
samples = rand(random_var_tuple,
                2, alg = RejectionSample)
# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,samples[2][N+x])
    push!(pred_vals,samples[2][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Simulation",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Pred"],
        lw = 1.25)

## Sample with conditional
samples = rand(random_var_tuple, random_var_tuple[7500] < 55,
                50, alg = RejectionSample)

# Histogram to test that conditional worked
compile = []
for x in 1:50
    push!(compile, samples[x][7500])
end
histogram(compile, bins = 10,
        title = "Prey at T=2500 (<55 Intervention)")

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,samples[1][N+x])
    push!(pred_vals,samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Prey <55 at T=2500 Simulation",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Pred"],
        lw = 1.25)

## Counterfactual
samples = rand(random_var_tuple, 1, alg = RejectionSample)
# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,samples[1][N+x])
    push!(pred_vals,samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Simulation",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Pred"],
        lw = 1.25)

new_prey_5 = replace(prey_list[5], prey_list[4] => 80.0)
deleteat!(prey_list, 5)
push!(prey_list, new_prey_5)
random_var_tuple = (Tuple(x for x in hazards_list)...,
                Tuple(x for x in prey_list)...,
                Tuple(x for x in pred_list)...,
                Tuple(Any[spawn_prey,prey2pred,pred_dies,theta])...)
samples = rand(new_prey_5, 1, alg = RejectionSample)

diffsamples = rand(random_var_tuple[7500],10)

## Test
x = normal(100, 10)
function deterministic(rng)
    if x(rng) > 100
        x(rng) + 10
    else
        x(rng) - 50
    end
end
deterministic_ = ciid(deterministic)
function deterministic_two(rng)
    deterministic_(rng) - 10
end
deterministic_two_ = ciid(deterministic_two)
tuple_test = (x, deterministic_, deterministic_two_)
samples = rand(tuple_test, tuple_test[3]<110, 100, alg = RejectionSample)
compile = []
for x in 1:100
    push!(compile, samples[x][3])
end
histogram(compile, bins = 10)

## Broken hazard
function test_hazards(rng)
    ecology = Dict("prey" => prey_init(rng), "pred" => pred_init(rng))

    hazards = Dict(
        "spawn_prey" => theta(rng)["spawn_prey"] * ecology["prey"],
        "prey2pred" => theta(rng)["prey2pred"] * ecology["prey"] * ecology["pred"],
        "pred_dies" => theta(rng)["pred_dies"] * ecology["pred"]
        )

    vals = collect(values(hazards))
    sum_vals = sum(vals)
    prob_vals = vals/sum_vals
    rand(Multinomial(1, prob_vals))[1]
end

hazards_temp = ciid(test_hazards)
rand((prey_init, pred_init, spawn_prey, prey2pred, pred_dies, theta,hazards_temp),
        hazards_temp<2, 100, alg = RejectionSample)
