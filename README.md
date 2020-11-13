# Causal_Inference_Project

## Installing Omega

I had some issues installing Omega (It wouldnt load after installation). I'd recommend installing it in the following sequence:

Pkg.add(Pkg.PackageSpec(;name="Flux", version="0.9.0"))

Pkg.add("Omega")

- Devon

## Literature

A Language for Counterfactual Generative Models - Overview of the Omega programming language recommended by Robert (most important)***

Stochastic.Modelling.for.Systems.Biology.Second.Edition - Info about Gillespie on Page 182 - 188, example of an SBML on pg 188 (more important)**

https://en.wikipedia.org/wiki/Gillespie_algorithm - Wiki for Gillespie is actually pretty helpful as well (more important)**


Leveraging Structured Biological Knowledge for Counterfactual Inference: a Case Study of Viral Pathogenesis - Paper written by previous students and Robert Ness.

Beyond Structural Causal Models Causal Constraints Models - Proposes CCMs as an alternative to SCMS when modeling dynamic systems (Less important)

Systems Biology: Basic Pathway Modeling Techniques - Intro to biochem kinetics and tellurium (Less important)

Generating Synthetic Signaling Networks - A paper that goes over how to generate random signaling networks (Less important)
