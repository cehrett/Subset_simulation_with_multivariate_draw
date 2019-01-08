# Subset simulation with multivariate draw
Carl Ehrett, 2017

## Nature of project
Subset simulation (SS) is a method of estimating low probability events. Suppose that one has a computer model predicting outcomes in a situation of interest. An example would be a model that predicts the amount of downtime suffered by a power plant following an earthquake. To assess the robustness of one's power infrastructure, one might wish to know the probability that the plant's downtime will exceed (e.g.) two months. Downtimes this long are extremely rare, and as a result, it is difficult to estimate their probability. 

For example, suppose that the probability of suffering such a long downtime following an earthquake (given that an earthquake of a given magnitude occurs) is (unbeknownst to us) 0.001. Thus our region of interest F is the set of downtimes exceeding two months, and we wish to estimate the probability of F. A common method of estimating the probability of a subregion of a computer model's range is to perform Monte Carlo estimation of that subregion, in which inputs to the model are randomly generated (from a prior on those inputs) and then the computer model output corresponding to those outputs is observed. Where M sets of inputs are drawn, some fraction p/M of the computer model outputs will be seen to fall in the region of interest F. Thus p/M becomes our estimate of the probability of F. This approach is problematic if our computer model is computationally expensive, since we will have to run the model on average 1000 times to get a single observation falling in the region F. So unless we are willing to spend copious time and computational resources on running our model millions of times, our estimate will have very high variance.

By contrast, subset simulation uses a series of intermediate failure regions, finding their conditional probabilities and using these to estimate the probability of F. For example, we might first find the probability of downtime exceeding one day, then the probability that downtime exceeds one week *given* that it exceeds one day, then the probability that downtime exceeds one month *given* that it exceeds one week, and so on. Markov chain Monte Carlo methods (MCMC) are used to sample from each of these conditional distributions. The product of the probabilities of these failure regions gives the probability of F, the region of interest.

Here I adapt SS to perform well with inputs that are correlated with respect to the failure region. That is, SS typically draws the model inputs when performing each MCMC using independent proposal densities. When the failure region F is spherical, this works well. But when F is not spherical, the inputs are correlated with respect to the failure region, and we can use the correlation observed in the intermediate failure regions to improve the proposal density of the model inputs and thereby massively improve the efficiency of the algorithm.

## Files
[demo.m] is a MATLAB script which demonstrates the proposed methodology and "stock" SS, using the example computer model tpf.m.

[tpf.m] is an example performance function to demonstrate SS. What tpf does is to find the inverse distance of a point from the origin, scaled using a hyperellipse. Different dimensionalities of inputs can be used, but by default demo.m uses a 100-dimensional input space of which four dimensions form the hyperellipse. So the set of points such that tpf returns a value greater than (e.g.) 1 is the set of points inside a shape which is a hyperellipse in four dimensions and which fills the unit hypercube in the remaining 96 dimensions. The only important point here is that having a tpf score above 1 means that you are inside a certain subregion of the unit hypercube, and the volume of that region is difficult to estimate via Monte Carlo because it is so small. Subset simulation can be used to estimate its volume.

[tpfparams.mat] contains settings used by the example performance function tpf.m.

[rot.m] is a script that produces a random rotation matrix of any desired dimension; it was used to produce tpfparams.mat.

[SS.m] performs subset simulation. To do so, it also calls MAS.m.

[MAS.m] performs the Metropolis-Hastings sampling for the MCMC in each intermediate failure region.

## To use this software

Load the file demo.m to run the demo. If desired, you can there replace tpf.m with whichever performance function you wish to use. the variable B in demo.m should also be changed to define whatever is the region of interest in the computer model you use to replace tpf.

Running the demo will show that the proposed version of SS can be significantly faster and more computationally efficient than the "stock" version.
