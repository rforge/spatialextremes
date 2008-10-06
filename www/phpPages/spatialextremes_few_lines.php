In this section, we explicit some of the most useful functions of
the package. However, for a full description, users may want to
have a look to the package vignette and the html help of the
package.

<h4 id="GPDComp">Fit a max-stable process - assuming unit Frechet
  margins:</h4>
<div class="Rcodes">
  <span class="Rcomments"> 
    ##Fit the the Smith's model to data:<br/>
  </span>
  <a href="#Smith" class="Routs">
    <code> fitmaxstab(data, coord, "gauss", fit.marge = FALSE) </code>
    <span class="Rcomments">
      <?php $tab = include "Routs/smith.out"; ?>
    </span>
  </a>
  <br/>
  <span class="Rcomments">
    ##Fit the Schlather's model with the powered exponential covariance
    function:<br/>
  </span>
  <a href="#Powexp" class="Routs">
    <code>fitmaxstab(data, corrd, "powexp", fit.marge = FALSE)</code>
    <span class="Rcomments">
      <?php $tab = include "Routs/powexp.out"; ?>
    </span>
  </a>
</div>
<h4 id="GPDFit">Fit a max-stable process - allowing for common GEV
  margins</h4>
<div class="Rcodes">
  <span class="Rcomments">
    ##With linear models for the GEV parameters:<br/>
  </span>
  <a  href="#linearModGEV" class="Routs">
    <code>fitmaxstab(data, coord, loc.form = y ~ lon, scale.form = y ~
      lat, shape.form = y ~ 1)</code>
  </a>
  <br/>
  <span class="Rcomments">
    ##With a p-spline on the location parameter:<br/> 
  </span>
  <a href="#splineLoc" class="Routs">
    <code>
      loc.form &lt;- y ~ rb(lon, knots, degree, penalty)<br/>
      fitmaxstab(data, coord, "gauss", loc.form = loc.form, scale.form =
      y ~ lat, shape.form = y ~ 1)</code>
  </a>
</div>  
<h4 id="GPDPlot">Plots</h4>
<div class="Rcodes">
  <span class="Rcomments">
    ##Covariance functions:<br/></span>
  <a href="images/covariance.png" class="Routs">
    <code class="Routs">covariance(sill = 1, range = 1, smooth = 0.5, cov.mod = "whitmat");</code>
    <span>
      <img src="images/covariance.png"
	   alt="covariance function"
	   height="300"/>
    </span>
  </a>
  <br/>
  <span class="Rcomments">
    ##Extremal coefficient:<br/> </span>
  <a class="Routs" href="images/extcoeff.png">
    <code class="Routs">extcoeff(fitted.object);</code>
    <span>
      <img src="images/extcoeff.png" 
	   alt="extremal coefficient"
	   height="300"/>
    </span>
  </a>
  <br/>
  <span class="Rcomments">
    ##Prediction of return levels:<br/> 
  </span>
  <a class="Routs" href="images/map.png">
    <code class="Routs">map(fitted.object, "quant", ret.per =
    50)</code>
    <span>
      <img src="images/map.png"
	   alt="prediction of return levels"
	   height="150"/>
    </span>
  </a>
  <br/>
  <a class="Routs" href="images/condmap.png">
    <code class="Routs">condmap(schlather, c(1, 1), seq(0, 10, length
    = 20), seq(0,10, length = 20))</code>
    <span>
      <img src="images/condmap.png"
	   alt="prediction of conditional return levels"
	   height="150"/>
    </span>
  </a>
  <br/>
  <span class="Rcomments">
    ##Semi-parametric estimates of the extremal coefficient:<br/> 
  </span>
  <a class="Routs" href="images/fitextcoeff.png">
    <code class="Routs">fitextcoeff(data, coord); fitextcoeff(data,
    coord, angles = seq(-pi, pi, length = 4))</code>
    <span>
      <img src="images/fitextcoeff.png"
	   alt="semi-parametric estimates of the extremal coefficient"
	   height="300"/>
    </span>
  </a>
  <br/>
  <span class="Rcomments">
    ##Profile composite likelihood:<br/> 
  </span>
  <a class="Routs" href="images/proflik.png">
    <code class="Routs">profile(fitted.object, "cov11", range = c(80,
    200))</code>
    <span>
      <img src="images/proflik.png"
	   alt="Profile composite likelihood"
	   height="300"/>
    </span>
  </a>
</div>
