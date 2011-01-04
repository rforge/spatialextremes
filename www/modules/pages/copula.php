<?php $nav_en_cours = 'learnmore'; $title = 'Copulas'?>
<?php include("modules/pages/menuLearnMore.php"); ?>
<div id="right">
  <div id="right_top"></div>
  <div id="right_bg">
    Before introducing the copula framework, I would like to clarify
    some points. I'm not "a big fan" of copulas especially when our
    interest is to model spatial extremes. Indeed the copula framework
    can be misleading since the spatial dependence of extremes might
    be falsely taken into account. The extreme value theory suggests
    that one should use max-stable copula and this corresponds
    actually to consider the finite dimensional distributions of a
    max-stable process. However I decided to implement copulas mainly
    for educational purposes.
    <h1>What are copulas?</h1>
    <div class="text">
      We introduce copulas by considering the most used copula: the
      Gaussian copula. Recall that we are interesting in modeling
      spatial extremes and in particular univariate arguments suggest
      that block maxima should be well described by a GEV
      distribution. If we
      denote <img src="http://latex.codecogs.com/gif.latex?F_x"> the
      distribution
      of <img src="http://latex.codecogs.com/gif.latex?Y(x)"> for
      all <img src="http://latex.codecogs.com/gif.latex?x \in
      \mathbb{R}^d"> one would write <br /><br />
      <div style="text-align: center">
	<img src="http://latex.codecogs.com/gif.latex?
		  \Pr[Y(x_1) \leq y_1, \leq, Y(x_k) \leq y_k] =
		  \Phi_\Sigma\{\Phi^{-1}(u_1), \ldots,
		  \Phi^{-1}(u_k)\},
		  ">
      </div>
      where <img src="http://latex.codecogs.com/gif.latex?\Phi_\Sigma">
      is the multivariate Normal distribution with zero mean and
      covariance
      matrix <img src="http://latex.codecogs.com/gif.latex?\Sigma">
      whose diagonal elements are all equal to
      unity, <img src="http://latex.codecogs.com/gif.latex?\Phi^{-1}">
      the quantile function of a standard Normal random variable
      and <img src="http://latex.codecogs.com/gif.latex?u_i =
      F_{x_i}(y_i)"> for
      all <img src="http://latex.codecogs.com/gif.latex? i = 1,
      \ldots, k">.

      Actually the above equation corresponds to the use of a Gaussian
      copula but other copulas can be used by taking for instance the
      multivariate Student distribution. Note that the Gaussian copula
      is asymptotically independent which implies that the extremes
      will occur independently from one location to another one ---
      which is not what we really want for spatial extreme don't we?
      The Student copula however is asymptotically dependent but tends
      to underestimate the spatial dependence of extreme events.
    </div>
    <h1>Function "fitcopula": Fit copulas</h1>
    <div class="text">
      The density corresponding to the Gaussian copula is easily found
      to be <br /><br />
      <div style="text-align: center">
	<img src="http://latex.codecogs.com/gif.latex?
		  f(y_1, \ldots, y_k) =
		  \frac{\varphi_\Sigma\{\Phi^{-1}(u_1), \ldots,
		  \Phi^{-1}(u_k)\}}{\prod_{i=1}^k
		  \varphi\{\Phi^{-1}(u_i)\}} \prod_{i=1}^k
		  f_{x_i}(y_i), ">
      </div>
      where <img src="http://latex.codecogs.com/gif.latex?\varphi_\Sigma">
      is the multivariate Normal density related to
      <img src="http://latex.codecogs.com/gif.latex?\Phi_\Sigma"> and
      <img src="http://latex.codecogs.com/gif.latex?f_{x_i}"> is the
      density related
      to <img src="http://latex.codecogs.com/gif.latex?F_{x_i}">. If
      one use the Student copula this would give a similar expression
      where the Normal densities and quantile functions are
      substituted for their Student's analogues.
    </div>
  </div>
  <div id="right_bot"></div>
</div>
