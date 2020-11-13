using Omega
using UnicodePlots

weight = β(2.0,2.0)
beta_samples = rand(weight, 10000)

UnicodePlots.histogram(beta_samples)

nflips = 4
coinflips_ = [bernoulli(weight, Bool) for i = 1:nflips]
