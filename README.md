# Causal Inference in Complex Systems Using Omega

This is an implementation of the Lotka-Volterra Predatar Prey model using Gillespie simulation with the probabilistic programming language Omega, allowing us to ask causal questions with the model

In this project we attempt to overcome some of the restrictions that come from using an SCM by applying probabilistic programming to the program
Probabilistic programming allows us to model the time and cycle component of systems, through the use of recursion and dynamically determined interventions (such as intervening when a value reaches some threshold)
In this project we use the probabilistic programming Omega, implemented in Julia
Once program is implemented, could scale it to any system that is described by a set of stochastic (or ordinary) differential equations

## Installing Omega

Pkg.add(Pkg.PackageSpec(;name="Flux", version="0.9.0"))

Pkg.add("Omega")


## Literature

A Language for Counterfactual Generative Models - Overview of the Omega programming language recommended by Robert (most important)***

Stochastic.Modelling.for.Systems.Biology.Second.Edition - Info about Gillespie on Page 182 - 188, example of an SBML on pg 188 (more important)**

https://en.wikipedia.org/wiki/Gillespie_algorithm - Wiki for Gillespie is actually pretty helpful as well (more important)**


Leveraging Structured Biological Knowledge for Counterfactual Inference: a Case Study of Viral Pathogenesis - Paper written by previous students and Robert Ness.

Beyond Structural Causal Models Causal Constraints Models - Proposes CCMs as an alternative to SCMS when modeling dynamic systems (Less important)

Systems Biology: Basic Pathway Modeling Techniques - Intro to biochem kinetics and tellurium (Less important)

Generating Synthetic Signaling Networks - A paper that goes over how to generate random signaling networks (Less important)
