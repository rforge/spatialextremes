<?php include("modules/pages/menuLearnMore.php"); ?>
<div id="right">
  <div id="right_top"></div>
  <div id="right_bg">
    <h1>What are Bayesian Hierarchical models?</h1>
    <div class="text">
      Hierarchical models have different layers of variations which
      must be modelled. When trying to model spatial extremes we can
      think of (at least) two layers: a layer that determines the
      marginal behaviour of extremes and another layer that controls
      the spatial dependence of extremes. Unfortunately because the
      likelihood of max-stable processes is not available it prevents
      the use of a spatial dependence layer. Consequently this type of
      models will never be able to model the dependence between say
      two weather stations. However its flexibility for modelling the
      marginal behaviour and will be of great interest if one is
      interested in modelling only the pointwise distribution of
      extreme events.

      <br /> If we dispose of block, e.g., annual, maxima observed at
      locations <img src="http://latex.codecogs.com/gif.latex?x_1,
      \ldots, x_k \in \mathbb{R}^d">, univariate extreme value
      arguments suggests that these block maxima might be modelled by
      a GEV distribution. The key idea of the latent process approach
      is to assume that the GEV parameters vary smoothly over space
      according to a stochastic
      process <img src="http://latex.codecogs.com/gif.latex?
      \{S(x)\}">. The SpatialExtremes package use Gaussian processes
      for this and assume that the Gaussian processes related to each
      GEV parameter are mutually independent. For instance we take
      <br /><br />
      <div style="text-align: center">
	<img src="http://latex.codecogs.com/gif.latex?\mu(x) = f_\mu(x;
		  \beta_\mu) + S_\mu(x; \psi_\mu),">
      </div>
      where <img src="http://latex.codecogs.com/gif.latex?f_\mu"> is a
      trend surface depending on regression parameters
      <img src="http://latex.codecogs.com/gif.latex?\beta_\mu">
      and <img src="http://latex.codecogs.com/gif.latex?S_\mu"> is a
      zero mean stationary Gaussian process with covariance function
      depending on
      parameters <img src="http://latex.codecogs.com/gif.latex?\psi_\mu">. Similar
      expressions might be used for the two remaining GEV parameters.
      <br /> Then conditional on the values of the three Gaussian
      processes at the weather stations, the block maxima are assumed
      to follow a GEV distribution <br /><br />
      <div style="text-align: center">
	<img src="http://latex.codecogs.com/gif.latex?Y_i(x_j) \mid
		  \{\mu(x_j); \sigma(x_j), \xi(x_) \} \sim \mbox{GEV}\{\mu(x_j),
		  \sigma(x_j), \xi(x_j)\}, \quad i = 1, \ldots, n, \quad j = 1,
		  \ldots, k,">
      </div>
      independently for each location.
    </div>
    <h1>Function "latent": Draw a Markov chain from the posterior
      distribution</h1>
    <div class="text">
      The full conditional distributions are easily shown to be
      <br /><br />
      <div style="text-align: center">
	<img src="http://latex.codecogs.com/gif.latex?
		  \pi(\mu \mid \cdots) &\propto& \pi(\mu \mid \beta_\mu,
		  \psi_\mu) \pi(y \mid \mu, \sigma, \xi)\\
		  \pi(\psi_\mu \mid \cdots) &\propto& \pi(\psi_\mu)
		  \pi(\mu \mid \beta_\mu, \psi)\\
		  \psi(\beta_\mu \mid \cdots) &\propto& \pi(\beta_\mu)
		  \pi(\mu \mid \beta_\mu, \psi_\mu),
		  ">
      </div>
      where <img src="http://latex.codecogs.com/gif.latex?
      \pi(\psi_\mu), \pi(\beta_\mu)"> are prior distributions. The
      full conditional distributions related to the two remaining GEV
      parameters give similar expressions.
    </div>
  </div>
  <div id="right_bot"></div>
</div>
