<?php
   $domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
   $group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']); 
   $themeroot='http://r-forge.r-project.org/themes/rforge/';
   echo '<?xml version="1.0" encoding="UTF-8"?>';
   ?>
<!DOCTYPE html
	  PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"> 
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en"> 
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8"/> 
    <title>The SpatialExtremes Package: An R package to model spatial extremes</title>
    <link rel="stylesheet" media="screen" type="text/css"
	  title="My Design" href="css/mystyle.css" />
    <link rel="shortcut icon" type="image/png" 
	  href="images/favicon.png" />
    <meta name="description" lang="en" 
	  content="The SpatialExtremes package aims to provide tools 
		   to model spatial extremes."/>
    <meta name="keywords" 
	  content="spatial extremes, extreme value theory, max-stable
		   processes, composite likelihood, latent variable,
		   copula"/>
    <script type="text/javascript"
    src="javascript/lightbox.js"></script>
    <script type="text/javascript">

      var _gaq = _gaq || [];
      _gaq.push(['_setAccount', 'UA-19518950-1']);
      _gaq.push(['_trackPageview']);
      
      (function() {
      var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
      ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
      var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
      })();
      
    </script>
  </head>
  <body>
    <div id="conteneur">
      <!-- header part -->
      <h1 id="header"><a href="#features" title="The SpatialExtremes Package">
	</a>
      </h1>
      <?php include "phpPages/menu.php"; ?>
      <!-- menu part --> 
      <!-- main matter part -->
      <div id="contenu">
	<?php include "phpPages/welcome.php"; ?>
	<h2 id="features">Features</h2>
	<?php include "phpPages/features.php"; ?>
	<h2 id="spatialextremes_few_lines">The SpatialExtremes Package
	  in a Few Lines</h2>
	<?php include "phpPages/spatialextremes_few_lines.php"; ?>
	<h2 id="manuals">Manuals</h2> 
	<?php include "phpPages/manuals.php"; ?>
	<h2 id="join">Contribute to the Project</h2> If you are
	interested in joining the project, you must first create an
	<a href="https://r-forge.r-project.org/account/register.php">R-forge
	  account</a>. Then, just
	go <a href="https://r-forge.r-project.org/tracker/admin/?group_id=265">
	  here</a>.
	<h2 id="contact">Contact</h2> Any suggestions, feature requests,
	bugs: <a href="https://r-forge.r-project.org/tracker/admin/?group_id=265">
	  select the appropriate tracker</a><br/> Author: Mathieu
	Ribatet <a href="http://people.epfl.ch/mathieu.ribatet">
	  (homepage)</a>
      </div>
    </div>
    <!-- footnote part -->
    <div id="footnote">
      <p>
	<a href="http://validator.w3.org/check?uri=referer">
	  <img src="http://www.w3.org/Icons/valid-xhtml10-blue"
	       alt="Valid XHTML 1.0 Transitional" height="31"
	       width="88"/>
	</a>
	<a href="http://validator.w3.org/check?uri=referer">
	  <img src="http://www.w3.org/Icons/valid-css"
	       alt="Valid CSS" height="31"
	       width="88"/>
	</a>
      </p>
    </div>
  </body>
</html>
