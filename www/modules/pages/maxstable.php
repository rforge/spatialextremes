<?php $nav_en_cours = 'learnmore'; $title = 'Max-stable processes'?>
<?php include("modules/pages/menuLearnMore.php"); ?>
<div id="right">
  <div id="right_top"></div>
  <div id="right_bg">
    <h1>What are max-stable processes?</h1>
    <div class="text">
    Max-stable processes are the extension of the multivariate extreme
    value theory to the infinite dimensional setting. More precisely
    consider a stochastic
    process <img src="http://latex.codecogs.com/gif.latex?\{Y(x)\}, x
    \in \mathbb{R}^d">, having continuous sample paths. Then the
    limiting process <br /><br />
    <div style="text-align: center">
      <img src="http://latex.codecogs.com/gif.latex?\left\{
		\max_{i=1, \ldots, n} \frac{Y_i(x) -b_n(x)}{a_n(x)}
		\right\}_{x \in \mathbb{R}^d} \longrightarrow \{Z(x)\}_{x
		\in \mathbb{R}^d}, \qquad n \to \infty,">
    </div>
    where <img src="http://latex.codecogs.com/gif.latex?Y_i"> are
    independent replications
    of <img src="http://latex.codecogs.com/gif.latex?Y">, <img src="http://latex.codecogs.com/gif.latex?a_n(x)>
    0"> and <img src="http://latex.codecogs.com/gif.latex?b_n(x) \in
    \mathbb{R}^d"> are sequences of continuous functions and the
    limiting process <img src="http://latex.codecogs.com/gif.latex?Z">
    is assumed to be non degenerate. de Haan [1994] shows that the
    class of the limiting processes
    <img src="http://latex.codecogs.com/gif.latex?Z(x)">
    corresponds to that of max-stable processes and hence emphasizes on
    their use to model spatial extremes.
    
    <br>
    Interestingly there are two different ways of charactaresing
    max-stable processes: these are known as the spectral
    characterisations. An especially useful special case of the
    characterisation of de Haan [1984], is <br /><br />
    <div style="text-align: center">
      <img src="http://latex.codecogs.com/gif.latex?
		Z(x) = \max_{i\geq 1} \xi_i f(x - U_i)">
    </div>
    where <img src="http://latex.codecogs.com/gif.latex? \{(\xi_i,
		    U_i)\}_{i \geq 1}"> are the points of a Poisson point
    process
    on <img src="http://latex.codecogs.com/gif.latex?
		 (0, \infty] \times \mathbb{R}^d"> with intensity
    measure <img src="http://latex.codecogs.com/gif.latex? 
		      \mbox{d$\Lambda$}(\xi,u) = \xi^{-2} \mbox{d$\xi$} \mbox{du}">,
    <img src="http://latex.codecogs.com/gif.latex? f"> is a
    probability density function
    on <img src="http://latex.codecogs.com/gif.latex?
    \mathbb{R}^d">. The process
    <img src="http://latex.codecogs.com/gif.latex? Z"> defined above
    is a max-stable process with unit Frechet
    margins. Taking <img src="http://latex.codecogs.com/gif.latex?
    f"> as the multivariate Normal density with zero mean and
    covariance matrix <img src="http://latex.codecogs.com/gif.latex?
    \Sigma"> gives the Smith model [Smith, 1990] for which the
    bivariate distribution function is <br /><br />
    <div style="text-align: center">
      <img src="http://latex.codecogs.com/gif.latex?
		\Pr[Z(x_1) \leq z_1, Z(x_2) \leq z_2] = \exp\left\{-
		\frac{1}{z_1} \Phi\left(\frac{a}{2} + \frac{1}{a} \log
		\frac{z_1}{z_2} \right) - \frac{1}{z_2}
		\Phi\left(\frac{a}{2} + \frac{1}{a} \log \frac{z_2}{z_1}
		\right) \right\},">
    </div>
    where <img src="http://latex.codecogs.com/gif.latex? \Phi"> is the
    standard normal distribution function
    and <img src="http://latex.codecogs.com/gif.latex? a^2 = (x_1 -
    x_2)^T \Sigma^{-1} (x_1 -
    x_2)">, <img src="http://latex.codecogs.com/gif.latex?  x_1, x_2
    \in \mathbb{R}^d.">
    
    <br /> Another very useful spectral characterisation for unit
    Frechet max-stable processes is [Schlather, 2002]<br /><br />
    <div style="text-align: center">
      <img src="http://latex.codecogs.com/gif.latex?
		Z(x) = \max_{i \geq 1} \xi_i Y_i(x),">
    </div>
    where <img src="http://latex.codecogs.com/gif.latex? \{\xi_i\}_{i
		    \geq 1}"> are the points of a Poisson point
		    process
		    on <img src="http://latex.codecogs.com/gif.latex?
		    (0, \infty]"> with intensity
		    measure <img src="http://latex.codecogs.com/gif.latex?
		    \mbox{d$\Lambda$}(\xi) = \xi^{-2} \mbox{d$\xi$}">
		    and
    <img src="http://latex.codecogs.com/gif.latex? Y_i"> are
    independent replications of a positive stochastic process having
    continuous sample
    paths <img src="http://latex.codecogs.com/gif.latex? Y"> such
    that <img src="http://latex.codecogs.com/gif.latex?
    \mathbb{E}[Y(x)] = 1"> for
    all <img src="http://latex.codecogs.com/gif.latex? x \in
    \mathbb{R}^d">.
    
    <br /> Currently there are several useful models based on
    Schlather's characterisation. The first model, the Schlather's
    model, is to use <img src="http://latex.codecogs.com/gif.latex?
    Y(x) = \sqrt{2 \pi} \max \{0, \varepsilon(x)\}">, where
    <img src="http://latex.codecogs.com/gif.latex? \varepsilon"> is a
    standard Gaussian process. This leads to the bivariate distribution
    function <br /><br />
    <div style="text-align: center">
      <img src="http://latex.codecogs.com/gif.latex?
		\Pr[Z(x_1) \leq z_1, Z(x_2) \leq z_2] =
		\exp\left[-\frac{1}{2} \left(\frac{1}{z_1} +
		\frac{1}{z_2} \right) \left(1 + \sqrt{1 - \frac{2 \{1 +
		\rho(h) \} z_1 z_2}{(z_1 +z_2)^2}} \right) \right]",>
    </div>
    
    Another possibility is to
    take <img src="http://latex.codecogs.com/gif.latex? Y(x) = \exp\{
    \sigma \varepsilon(x) - \sigma^2 / 2\}">,
    where <img src="http://latex.codecogs.com/gif.latex?  \sigma >
    0">. This is known as the Geometric gaussian model for which the
    bivariate distribution is similar to the Smith model
    with <img src="http://latex.codecogs.com/gif.latex? a^2 = 2
    \sigma2 \{1 - \rho(h)\}">. A last possibility, which is a
    generalisation of the geometric Gaussian model, is to
    take <img src="http://latex.codecogs.com/gif.latex? Y(x) =
    \exp\{\varepsilon(x) - \gamma(x)\}">, where
    <img src="http://latex.codecogs.com/gif.latex?  \varepsilon"> is a
	      zero mean Gaussian process having stationary increments
	      and (semi)variogram
    <img src="http://latex.codecogs.com/gif.latex?  \gamma"> such
	      that <img src="http://latex.codecogs.com/gif.latex?
	      \varepsilon(o) = 0"> almost surely. Its bivariate
	      distribution function is again similar to the Smith
	      model
	      with <img src="http://latex.codecogs.com/gif.latex?  a^2
	      = 2 \gamma(x_1 - x_2)">.
    </div>
    
    <h1>Function "fitmaxstab": Fit max-stable processes</h1>
    <div class="text">
    Because max-stable processes are an extension of the multivariate
    extreme value theory to the infinite dimensional case when trying
    to fit such processes to we are faced to the same problem has in
    the finite dimensional setting.  <br />
    
    Suppose we have observed the
    process <img src="http://latex.codecogs.com/gif.latex?  \{Z(x)\}">
    at <img src="http://latex.codecogs.com/gif.latex?  x_1, \ldots,
    x_k \in \mathbb{R}^d">, then the multivariate distribution is of
    the form <br /><br />
    <div style="text-align: center">
      <img src="http://latex.codecogs.com/gif.latex?
		\Pr[Z(x_1) \leq z_1, \ldots, Z(x_k) \leq z_k] =
		\exp\left\{-V(z_1, \ldots, z_k) \},">
    </div>
    where <img src="http://latex.codecogs.com/gif.latex? V"> is a
    function having some homogeneity property called the exponent
    measure. Consequently if one wants to use the maximum likelihood
    estimator this will yield to a combinatorial explosion and the
    likelihood will be intractable. Instead one can have resort to
    composite likelihood estimator and an especially convenient choice
    is to maximize the pairwise log-likelihood <br /><br />
    <div style="text-align: center">
      <img src="http://latex.codecogs.com/gif.latex?
		\ell_p(\psi) = \sum_{i=1}^n \sum_{1\leq j < l \leq k}
		w_{j,l} \log f(z_{i,j}, z_{i,l}; \psi),">
    </div>
    where <img src="http://latex.codecogs.com/gif.latex?
		    \{w_{j,l}\}_{j,l}"> is a set of appropriate
		    weights
		    and <img src="http://latex.codecogs.com/gif.latex?
		    f(z_{i,j}, z_{i,l}; \psi)"> is the bivariate
		    density of a max-stable process --- for which
		    closed forms exist.
    
    Under some regularity conditions and in particular
    if <img src="http://latex.codecogs.com/gif.latex?  \psi"> is
    identifiable from the bivariate densities, then <br /><br />
    <div style="text-align: center">
      <img src="http://latex.codecogs.com/gif.latex?
		\sqrt{n} \{H(\psi_0) J(\psi_0)^{-1} H(\psi_0)\}^{1/2}
		(\hat{\psi} - \psi_0) \stackrel{\rm d}{\longrightarrow}
		N(0, \mbox{Id}), \qquad n \to \infty,">
    </div>
    where <img src="http://latex.codecogs.com/gif.latex?
		    H(\psi_0) = \mathbb{E}[\nabla^2 \ell_p(\psi_0)]"> and 
    <img src="http://latex.codecogs.com/gif.latex?
	      J(\psi_0) = \mbox{Var}[\nabla \ell_p(\psi_0)]">.
    </div>
    <h1>Function "rmaxstab": Simulate max-stable processes</h1>
    <div class="text">
      Well we can simulate a max-stable process using its limiting
    characterisation i.e. find the normalizing
    functions <img src="http://latex.codecogs.com/gif.latex?a_n(x)">
    and <img src="http://latex.codecogs.com/gif.latex?b_n(x)"> and
    simulate <img src="http://latex.codecogs.com/gif.latex?n">
    indendant replications of a stochastic
    process <img src="http://latex.codecogs.com/gif.latex?\{Y(x)\}">
    for which the normalizing function are based. This naive
    simulation technique has the disadvantage of requiring a large
    number of independent replicates but also pose the question of how
    many replicates should be simulated? A better strategy consists in
    using one of the spectral representation. In most cases this will
    be more efficient since one can simulate the points of a Poisson
    point process with
    intensity <img src="http://latex.codecogs.com/gif.latex?
    \mbox{d$\Lambda$}(\xi) = \xi^{-2} \mbox{d$\xi$}"> can be generated
    as follows <br /><br />
    <div style="text-align: center">
      <img src="http://latex.codecogs.com/gif.latex?
		\xi_i = \frac{1}{\sum_{j=1}^i E_j}, \qquad E_j
		\stackrel{\rm iid}{\sim} \mbox{Exp}(1).">
    </div>
    Now as standard exponential random variables are positive,
    <img src="http://latex.codecogs.com/gif.latex? \xi_i \to 0"> as
    <img src="http://latex.codecogs.com/gif.latex? i \to \infty">. If we
    further assume that the stochastic process used in the spectral
    characterisation is bounded above, we will require only a finite
    number of independent replications. If the process is not bounded then
    Schlather [2002] shows that it is however possible to get accurate
    simulations. <br />

    The methodology introduced above was the first approach to
    generate (approximate) realizations from a max-stable process. It
    is now possible to get exact simulation using the methodology
    developped by <a
    href="https://academic.oup.com/biomet/article/103/2/303/1744000">C. Dombry,
    S. Engelke and M. Oesting</a>. The package uses such a strategy whenever it is possible (or implemented ;-) )
    
  </div>
  </div>
  <div id="right_bot"></div>
</div>
