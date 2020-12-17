# Causal Inference in Complex Systems Using Omega

This is an implementation of the Lotka-Volterra Predatar Prey model using Gillespie simulation with the probabilistic programming language Omega, allowing us to ask causal questions with the model

In this project we attempt to overcome some of the restrictions that come from using an SCM by applying probabilistic programming to the program
Probabilistic programming allows us to model the time and cycle component of systems, through the use of recursion and dynamically determined interventions (such as intervening when a value reaches some threshold)
In this project we use the probabilistic programming Omega, implemented in Julia
Once program is implemented, could scale it to any system that is described by a set of stochastic (or ordinary) differential equations

Implement causal modeling on dynamic systems using the probabilistic programming language Omega, allowing us to do inference and counterfactuals. We will do this by taking a set of ODEs/SDEs and simulating data using the Gillespie algorithm. To start we will implement the Lotka-Volterra predator-prey model. Ideally the function will be extended to intake any SBML and run simulations on the model. We need to implement Gillespie in a way that we can make interventions at any time during the simulation. Once Gillespie is implemented we can make inferences and counterfactuals on the individual simulations. We want to look into both inferring rates and species abundances. Omega allows us to move away from DAGs and model more complex systems with complexities such as cycles and time periods.

Additionally, as a long shot goal, we would like to make comparisons between using probabilistic programming and Omega to using an SCM or other non-causal modeling methods. We would like to explore the exact benefits of using probabilistic programming.

## Deliverables

1.       Implementation of Gillespie in Julia
    
2.       Simulation the Lotka-Volterra model using Gillespie

3.       Implementation of the the simulation in Omega

4.       Run different inferences/counterfactuals using Omega
     

## How to run the project


### Prerequisites 

Julia 

Installation:
For Windows:
https://julialang.org/downloads/platform/#windows
For macOS: 
https://julialang.org/downloads/platform/#macos
For Linux:
https://julialang.org/downloads/platform/#linux_and_freebsd


The following julia packages must be installed to run the Jupyter notebooks in this project:

```
Omega
StatsBase
Random
Plots
Distributions
```
If you have problems installing Omega try:

Pkg.add(Pkg.PackageSpec(;name="Flux", version="0.9.0"))
Pkg.add("Omega")



## Authors

* [**Devon Kohler**]()

* [**Vartika Tewari**]()

* [**Ritwik Ananad**]()

* [**Derek Grubis**]()



## License

This project is licensed under Apache License 2.0 - see the [LICENSE.md](https://github.com/devonjkohler/Causal_Inference_Project/blob/main/LICENSE.md) file for details



## Literature

A Language for Counterfactual Generative Models - Overview of the Omega programming language 

Stochastic.Modelling.for.Systems.Biology.Second.Edition - Info about Gillespie on Page 182 - 188, example of an SBML on pg 188 

https://en.wikipedia.org/wiki/Gillespie_algorithm 


Leveraging Structured Biological Knowledge for Counterfactual Inference: a Case Study of Viral Pathogenesis 

Beyond Structural Causal Models Causal Constraints Models - Proposes CCMs as an alternative to SCMS when modeling dynamic systems 

Systems Biology: Basic Pathway Modeling Techniques - Intro to biochem kinetics and tellurium 

Generating Synthetic Signaling Networks - A paper that goes over how to generate random signaling networks 
