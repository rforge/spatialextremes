<?php
/**
 * Piwik - Open source web analytics
 * 
 * @link http://piwik.org
 * @license http://www.gnu.org/licenses/gpl-3.0.html Gpl v3 or later
 * @version $Id: Controller.php 526 2008-06-25 23:57:04Z matt $
 * 
 * @package Piwik_SitesManager
 */


/**
 * 
 * @package Piwik_SitesManager
 */
class Piwik_SitesManager_Controller extends Piwik_Controller
{
	function index()
	{
		$view = new Piwik_View('SitesManager/templates/SitesManager.tpl');
		
		$sites = Piwik_SitesManager_API::getSitesWithAdminAccess();
		foreach($sites as &$site)
		{
			$site['alias_urls'] = Piwik_SitesManager_API::getSiteUrlsFromId($site['idsite']);
		}
//		var_dump($sites);exit;
		$view->sites = $sites;
		echo $view->render();
	}
	
	function displayJavascriptCode()
	{
		$jsTag = Piwik::getJavascriptCode(Piwik_Common::getRequestVar('idsite',1), Piwik_Url::getCurrentUrlWithoutFileName());

		$view = new Piwik_View('SitesManager/templates/DisplayJavascriptCode.tpl');
		$view->jsTag = $jsTag;
		
		echo $view->render();
	}
}