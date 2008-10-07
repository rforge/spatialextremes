<?php
/**
 * Piwik - Open source web analytics
 * 
 * @link http://piwik.org
 * @license http://www.gnu.org/licenses/gpl-3.0.html Gpl v3 or later
 * @version $Id: API.php 518 2008-06-09 00:44:10Z matt $
 * 
 * @package Piwik_UserSettings
 */


require_once "DataFiles/Browsers.php";
require_once "DataFiles/OS.php";
		
/**
 * 
 * @package Piwik_UserSettings
 */
class Piwik_UserSettings_API extends Piwik_Apiable
{
	static private $instance = null;
	static public function getInstance()
	{
		if (self::$instance == null)
		{            
			$c = __CLASS__;
			self::$instance = new $c();
		}
		return self::$instance;
	}
	
	public function getResolution( $idSite, $period, $date )
	{
		Piwik::checkUserHasViewAccess( $idSite );
		$archive = Piwik_Archive::build($idSite, $period, $date );
		$dataTable = $archive->getDataTable('UserSettings_resolution');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ReplaceColumnNames');
		return $dataTable;
	}

	public function getConfiguration( $idSite, $period, $date )
	{
		Piwik::checkUserHasViewAccess( $idSite );
		
		$archive = Piwik_Archive::build($idSite, $period, $date );
		$dataTable = $archive->getDataTable('UserSettings_configuration');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackReplace', array('label', 'Piwik_getConfigurationLabel'));
		$dataTable->queueFilter('Piwik_DataTable_Filter_ReplaceColumnNames');
		return $dataTable;
	}

	public function getOS( $idSite, $period, $date )
	{
		Piwik::checkUserHasViewAccess( $idSite );
		$archive = Piwik_Archive::build($idSite, $period, $date );
		$dataTable = $archive->getDataTable('UserSettings_os');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ReplaceColumnNames');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackAddMetadata', array('label', 'logo', 'Piwik_getOSLogo'));
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackAddMetadata', array( 'label', 'shortLabel', 'Piwik_getOSShortLabel') );
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackReplace', array( 'label', 'Piwik_getOSLabel') );
		return $dataTable;
	}
		
	public function getBrowser( $idSite, $period, $date )
	{
		Piwik::checkUserHasViewAccess( $idSite );
		$archive = Piwik_Archive::build($idSite, $period, $date );
		$dataTable = $archive->getDataTable('UserSettings_browser');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ReplaceColumnNames');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackAddMetadata', array('label', 'logo', 'Piwik_getBrowsersLogo'));
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackAddMetadata', array('label', 'shortLabel', 'Piwik_getBrowserShortLabel'));
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackReplace', array('label', 'Piwik_getBrowserLabel'));
		return $dataTable;
	}
	
	public function getBrowserType( $idSite, $period, $date )
	{
		Piwik::checkUserHasViewAccess( $idSite );
		$archive = Piwik_Archive::build($idSite, $period, $date );
		$dataTable = $archive->getDataTable('UserSettings_browserType');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ReplaceColumnNames');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackAddMetadata', array('label', 'shortLabel', 'ucfirst'));
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackReplace', array('label', 'Piwik_getBrowserTypeLabel'));
		return $dataTable;
	}
	
	public function getWideScreen( $idSite, $period, $date )
	{
		Piwik::checkUserHasViewAccess( $idSite );
		$archive = Piwik_Archive::build($idSite, $period, $date );
		$dataTable = $archive->getDataTable('UserSettings_wideScreen');	
		$dataTable->queueFilter('Piwik_DataTable_Filter_ReplaceColumnNames');		
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackAddMetadata', array('label', 'logo', 'Piwik_getScreensLogo'));
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackReplace', array('label', 'ucfirst'));
		return $dataTable;
	}
	
	public function getPlugin( $idSite, $period, $date )
	{
		Piwik::checkUserHasViewAccess( $idSite );
		$archive = Piwik_Archive::build($idSite, $period, $date );
		$dataTable = $archive->getDataTable('UserSettings_plugin');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ReplaceColumnNames');		
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackAddMetadata', array('label', 'logo', 'Piwik_getPluginsLogo'));
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackReplace', array('label', 'ucfirst'));
		return $dataTable;
	}	
}

function Piwik_getPluginsLogo( $oldLabel )
{
	return  "plugins/UserSettings/images/plugins/". $oldLabel . ".gif";
}

function Piwik_getOSLabel($oldLabel)
{
	if(isset($GLOBALS['Piwik_Oslist_IdToLabel'][$oldLabel]))
	{
		return $GLOBALS['Piwik_Oslist_IdToLabel'][$oldLabel];
	}
	return 'UNK';
}


function Piwik_getOSShortLabel($oldLabel)
{
	if(isset($GLOBALS['Piwik_Oslist_IdToShortLabel'][$oldLabel]))
	{
		return $GLOBALS['Piwik_Oslist_IdToShortLabel'][$oldLabel];
	}
	return 'UNK';
}

function Piwik_getBrowserTypeLabel($oldLabel)
{
	if(isset(Piwik_UserSettings::$browserType_display[$oldLabel]))
	{
		return Piwik_UserSettings::$browserType_display[$oldLabel];
	}
	return Piwik_Translate('General_Unknown');
}


function Piwik_getConfigurationLabel($str)
{
	$values = explode(";", $str);
	
	$os = Piwik_getOSLabel($values[0]);
	$name = $values[1];
	$browser = 'Unknown';
	if(isset($GLOBALS['Piwik_BrowserList_IdToLabel'][$name]))
	{
		$browser = $GLOBALS['Piwik_BrowserList_IdToLabel'][$name];
	}
	
	$resolution = $values[2];
	
	return $os . " / " . $browser . " / " . $resolution;
}

function Piwik_getBrowserLabel($oldLabel)
{
	$name = Piwik_getBrowserId($oldLabel);
	$version = Piwik_getBrowserVersion($oldLabel);
	if(isset($GLOBALS['Piwik_BrowserList_IdToLabel'][$name]))
	{
		return $GLOBALS['Piwik_BrowserList_IdToLabel'][$name] . " ". $version;
	}
	return 'UNK';
}

function Piwik_getBrowserShortLabel($oldLabel)
{
	$name = Piwik_getBrowserId($oldLabel);
	$version = Piwik_getBrowserVersion($oldLabel);
	if(isset($GLOBALS['Piwik_BrowserList_IdToShortLabel'][$name]))
	{
		return $GLOBALS['Piwik_BrowserList_IdToShortLabel'][$name] . " ". $version;
	}
	return 'UNK';
}

function Piwik_getBrowserId($str)
{
	return substr($str, 0, strpos($str, ';'));
}

function Piwik_getBrowserVersion($str)
{
	return substr($str, strpos($str, ';') + 1);
}

function Piwik_getBrowsersLogo($label)
{
	$id = Piwik_getBrowserId($label);
	return  "plugins/UserSettings/images/browsers/". $id . ".gif";
}

function Piwik_getOSLogo($label)
{
	$path = "plugins/UserSettings/images/os/". $label . ".gif";
	return $path;
}

function Piwik_getScreensLogo($label)
{
	return "plugins/UserSettings/images/screens/" . $label . ".gif";
}