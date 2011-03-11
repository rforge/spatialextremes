<?php $nav_en_cours = 'learnmore';
      $title = 'Conditional simulation of max-stable processes'?>
<?php include("modules/pages/menuLearnMore.php"); ?>
<div id="right">
  <div id="right_top"></div>
  <div id="right_bg">
    <h1>What is a max-linear model?</h1>
    <div class="text">
      A unit Frechet max-linear process is defined by
      <br /><br />
      <div style="text-align: center">
	<img src="http://latex.codecogs.com/gif.latex?Z(x) =
      \max_{j=1, \ldots, p} f_j(x) Z_j, \qquad x \in \mathbb{R}^d,">
      </div>
      where <img src="http://latex.codecogs.com/gif.latex?Z_j"> are
      independent unit Frechet random variables
      and <img src="http://latex.codecogs.com/gif.latex?f_j"> are
      non-negative deterministic functions such that
      <br /><br />
      <div style="text-align: center">
	<img src="http://latex.codecogs.com/gif.latex?\sum_{j=1}^p
      f_j(x) = 1, \qquad x \in \mathbb{R}^d.">
      </div>

      <br /> Any unit Frechet max-stable process can be arbitrarly
      well approximated by a max-linear model by taking sufficiently
      large <img src="http://latex.codecogs.com/gif.latex?p"> and
      suitable
      functions <img src="http://latex.codecogs.com/gif.latex?f_j">.
    </div>
    <h1>Function "condrmaxstab": Performs conditional simulations for
      max-stable processes</h1>
    <div class="text">
      Since the
      functions <img src="http://latex.codecogs.com/gif.latex?f_j">
      are deterministics, Wang and Stoev [2011] proposed an efficient
      algorithm to generate conditional simulations from a max-linear
      model and thus approximate conditional simulations from a
      max-stable process.

      <br />However this approaches has some drawbacks since it is not
      clear how to find appropriate
      functions <img src="http://latex.codecogs.com/gif.latex?f_j">. Therefore
      the only available model is currently the (discretized) Smith
      model for which
      <br /><br />
      <div style="text-align: center">
	<img src="http://latex.codecogs.com/gif.latex?f_j(x) =
      \varphi(x - u_j ; \Sigma), \qquad x \in \mathbb{R}^d,">
      </div>
      where <img src="http://latex.codecogs.com/gif.latex?\varphi(\cdot;\Sigma)">
      is the zero mean (multivariate) normal density with covariance
      matrix <img src="http://latex.codecogs.com/gif.latex?\Sigma">
      and <img src="http://latex.codecogs.com/gif.latex?u_j"> are
      appropriately chose points
      in <img src="http://latex.codecogs.com/gif.latex?\mathbb{R}^d">.
    </div>
  </div>
  <div id="right_bot"></div>
</div>
