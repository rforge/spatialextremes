<!-- header begins -->
<div id="header">
  <div id="logo">
    <a href= "http://spatialextremes.r-forge.r-project.org">
      <div style= "width:290px; height:20px;
		   position:absolute;left:15%;top:75px;">
      </div>
    </a>
    <h2><a href="http://www.r-project.org/" id="R">An
      	R package for spatial extreme modelling</a></h2>
  </div>
  <div id="buttons">
    <a href="index.php" class="but"  title="">Home</a>
    <a href="index.php?module=pages&amp;action=learnmore"
       <?php if ($nav_en_cours == 'learnmore') {echo 'id="current"';}
	     else {echo 'class="but"';} ?> title="">Learn More</a>
    <a href="index.php?module=pages&amp;action=gallery"
       <?php if ($nav_en_cours == 'gallery') {echo 'id="current"';}
	     else {echo 'class="but"';} ?>
       title="">Gallery</a>
    <a href="index.php?module=pages&amp;action=aboutMe"
       <?php if ($nav_en_cours == 'aboutMe') {echo 'id="current"';}
	     else {echo 'class="but"';} ?>
    title="">Contact</a>
  </div>
</div>
<!-- header ends -->
