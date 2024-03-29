## Load Packages
using Omega
using StatsBase
using Random
using Plots
using Distributions


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

# Initiate starting values randomely
prey_init = normal(10., .001)
pred_init = normal(10., .001)

# Random variables for rates
spawn_prey = normal(1.5, .01)
prey2pred = normal(1.0, .0001)
pred_dies = normal(3.0, .0075)

# prey_init = normal(50., .001)
# pred_init = normal(100., .001)
#
# # Random variables for rates
# spawn_prey = normal(.9, .01)
# prey2pred = normal(.004, .0001)
# pred_dies = normal(.4, .0075)

## Generate Model
hazards_list = Any[]
prey_list = Any[]
pred_list = Any[]
Random.seed!(1234)
theta = ciid(generate_rates)
push!(prey_list, prey_init)
push!(pred_list, pred_init)
N = 1000

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
Random.seed!(1234)
samples = rand(random_var_tuple,
                1, alg = RejectionSample)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,samples[1][N+x])
    push!(pred_vals,samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Omega Simulation",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Predators"],
        lw = 1.25)

a = shuffle(prey_vals)
b = sortperm(a)
condition_index = b[490:499]

## Sample with conditional
hazards_list = Any[]
prey_list = Any[]
pred_list = Any[]
theta = ciid(generate_rates)
push!(prey_list, prey_init)
push!(pred_list, pred_init)
N = 500
for f in 2:N
    last = f - 1
    hazards_temp = ciid(get_hazards, last) # individual step
    if f in 150:151
        prey_temp_ = ciid(one_simulation_prey, last, transitions) # individual step
        prey_temp = cond(prey_temp_, prey_temp_ ==ₛ 5)
    else
        prey_temp = ciid(one_simulation_prey, last, transitions) # individual step
    end
    #if f in 500:750
    #    pred_temp_ = ciid(one_simulation_pred, last, transitions) # individual step
    #    pred_temp = cond(pred_temp_, pred_temp_ < 120)
    #else
    pred_temp = ciid(one_simulation_pred, last, transitions) # individual step
    #end
    push!(hazards_list, hazards_temp)
    push!(prey_list, prey_temp)
    push!(pred_list, pred_temp)
end

random_var_tuple = (Tuple(x for x in hazards_list)...,
                Tuple(x for x in prey_list)...,
                Tuple(x for x in pred_list)...,
                Tuple(Any[spawn_prey,prey2pred,pred_dies,theta])...)

Random.seed!(1234)
samples = rand(random_var_tuple,
                1, alg = SSMH)

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

# Histogram to test that conditional worked
compile = []
for x in 1:50
    push!(compile, samples[x][7500])
end
histogram(compile, bins = 10,
        title = "Prey at T=2500 (<55 Conditional)")
## Counterfactual
Random.seed!(1234)
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

replace_index = 16
new_prey = replace(prey_list[replace_index],
                    prey_list[replace_index - 1] => 10.0)
prey_list_replace = prey_list
prey_list_replace[replace_index] = new_prey

random_var_tuple_int = (Tuple(x for x in hazards_list)...,
                Tuple(x for x in prey_list_replace)...,
                Tuple(x for x in pred_list)...,
                Tuple(Any[spawn_prey,prey2pred,pred_dies,theta])...)
Random.seed!(1234)
adj_samples = rand(random_var_tuple_int, 1, alg = RejectionSample)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,adj_samples[1][N+x])
    push!(pred_vals,adj_samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Intervention Initial Prey=150",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Predators"],
        lw = 1.25)

## Int fix test
replace_index = 16
test = replace(prey_list[replace_index],
                    prey_list[replace_index - 1] => normal(30,1))

Random.seed!(1234)
samples = rand(random_var_tuple, 1, alg = RejectionSample)
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,samples[1][N+x])
    push!(pred_vals,samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Intervention Prey T10=300",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Pred"],
        lw = 1.25)

## Intervene on initial value
prey_init_new = replace(prey_init, prey_init => 150.)
prey_list_replace = prey_list
prey_list_replace[1] = prey_init_new
random_var_tuple_int = (Tuple(x for x in hazards_list)...,
                Tuple(x for x in prey_list_replace)...,
                Tuple(x for x in pred_list)...,
                Tuple(Any[spawn_prey,prey2pred,pred_dies,theta])...)
Random.seed!(1234)
adj_samples = rand(random_var_tuple_int, 1, alg = RejectionSample)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,adj_samples[1][N+x])
    push!(pred_vals,adj_samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Intervention Initial Prey=150",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Predators"],
        lw = 1.25)

## Recreate Plots
hazards_list = Any[]
prey_list = Any[]
pred_list = Any[]
Random.seed!(2)
theta = ciid(generate_rates)
push!(prey_list, prey_init)
push!(pred_list, pred_init)
N = 200

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

## A
Random.seed!(22)
samples = rand(random_var_tuple,
                1, alg = RejectionSample)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,samples[1][N+x])
    push!(pred_vals,samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Omega Simulation",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Predators"],
        lw = 1.25)

## B
condition_index = [27, 75, 122, 187] #27

hazards_list = Any[]
prey_list = Any[]
pred_list = Any[]
Random.seed!(2)
theta = ciid(generate_rates)
push!(prey_list, prey_init)
push!(pred_list, pred_init)
N = 200
for f in 2:N
    last = f - 1
    hazards_temp = ciid(get_hazards, last) # individual step
    if f in condition_index
        prey_temp_ = ciid(one_simulation_prey, last, transitions) # individual step
        prey_temp = cond(prey_temp_, prey_temp_ >= 6.)
    else
        prey_temp = ciid(one_simulation_prey, last, transitions) # individual step
    end
    pred_temp = ciid(one_simulation_pred, last, transitions) # individual step
    push!(hazards_list, hazards_temp)
    push!(prey_list, prey_temp)
    push!(pred_list, pred_temp)
end

random_var_tuple = (Tuple(x for x in hazards_list)...,
                Tuple(x for x in prey_list)...,
                Tuple(x for x in pred_list)...,
                Tuple(Any[spawn_prey,prey2pred,pred_dies,theta])...)

Random.seed!(22)
samples = rand(random_var_tuple,
                1, alg = RejectionSample)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals, samples[1][N+x])
    push!(pred_vals, samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Condition On Too Many Prey",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Pred"],
        lw = 1.25)

## C/E
replace_index = 15
new_prey = replace(prey_list[replace_index],
                    prey_list[replace_index - 1] => 10.0)
prey_list_replace = prey_list
prey_list_replace[replace_index] = new_prey
random_var_tuple_int = (Tuple(x for x in hazards_list)...,
                Tuple(x for x in prey_list_replace)...,
                Tuple(x for x in pred_list)...,#pred_list_replace
                Tuple(Any[spawn_prey,prey2pred,pred_dies,theta])...)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,adj_samples[1][N+x])
    push!(pred_vals,adj_samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Pred=6 at t=15",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Predators"],
        lw = 1.25)

function increase_prey()

    Random.seed!(rand(1:10000))
    samples = rand(random_var_tuple, 1, alg = RejectionSample)
    adj_samples = rand(random_var_tuple_int, 1, alg = RejectionSample)

    # extract run results and plot
    prey_norm_vals = []
    prey_adj_vals = []
    for x in 1:(N-1)
        push!(prey_norm_vals,samples[1][N+x])
        push!(prey_adj_vals,adj_samples[1][(N*2)+x])
    end

    return sum(prey_adj_vals) - sum(prey_norm_vals)
end

check = [increase_prey() for x=1:100]

histogram(check, bins = 15,
        title = "blah")

##D/F
replace_index = 15
new_pred = replace(pred_list[replace_index],
                    pred_list[replace_index - 1] => 12.0)
pred_list_replace = pred_list
pred_list_replace[replace_index] = new_pred
random_var_tuple_int = (Tuple(x for x in hazards_list)...,
                Tuple(x for x in prey_list)...,
                Tuple(x for x in pred_list_replace)...,#pred_list_replace
                Tuple(Any[spawn_prey,prey2pred,pred_dies,theta])...)

Random.seed!(22)
adj_samples = rand(random_var_tuple_int, 1, alg = RejectionSample)

# extract run results and plot
prey_vals = []
pred_vals = []
for x in 1:(N-1)
    push!(prey_vals,adj_samples[1][N+x])
    push!(pred_vals,adj_samples[1][(N*2)+x])
end

plot(hcat(prey_vals,pred_vals),
        title = "Pred=6 at t=15",
        xlabel = "Time",
        ylabel = "Quantity",
        label = ["Prey" "Predators"],
        lw = 1.25)


function increase_pred()

    samples = rand(random_var_tuple, 1, alg = RejectionSample)
    adj_samples = rand(random_var_tuple_int, 1, alg = RejectionSample)

    # extract run results and plot
    prey_norm_vals = []
    prey_adj_vals = []
    for x in 1:(N-1)
        push!(prey_norm_vals,samples[1][N+x])
        push!(prey_adj_vals,adj_samples[1][(N*2)+x])
    end

    return sum(prey_adj_vals) - sum(prey_norm_vals)
end

check = [increase_pred() for x=1:10]

histogram(check,
        title = "blah")
