<?php
   $domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
   $group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']); 
   $themeroot='http://r-forge.r-project.org/themes/rforge/';
   ?>

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="content-type" content="text/html; charset=utf-8" />
    <title>SpatialExtremes: An R package for modelling spatial
    extremes</title>
    <meta name="keywords" content="spatial extremes, extreme value
    theory, max-stable processes, composite likelihood, latent
    variable, copula, R package, Bayesian hierarchical models,
    geostatistics, extremes" />
    <meta name="description" content="The SpatialExtremes package aims
    to provide tools  to model spatial extremes." />
    <meta name="Author" content="Mathieu Ribatet" />
    <link href="/styles.css" rel="stylesheet" type="text/css" />
    <link rel="stylesheet" href="imageflow.css" type="text/css" />
    <link rel="stylesheet" href="highslide.css" type="text/css" />
    <script src="js/jquery-1.4.2.min.js" type="text/javascript"></script>
    <script src="js/imageflow.js" type="text/javascript"></script>
    <script src="js/highslide.js" type="text/javascript"></script>
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
    <div id="foot_bg">
      <div id="main">
	<?php include 'global/header.php'; ?>
	<!-- content begins -->
       	<div id="content">
