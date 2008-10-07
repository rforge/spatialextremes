<?php
/**
 * Piwik - Open source web analytics
 * 
 * @link http://piwik.org
 * @license http://www.gnu.org/licenses/gpl-3.0.html Gpl v3 or later
 * @version $Id: API.php 561 2008-07-21 00:00:35Z matt $
 * 
 * @package Piwik_VisitorInterest
 */


/**
 * 
 * @package Piwik_VisitorInterest
 */
class Piwik_VisitorInterest_API extends Piwik_Apiable
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

	public function getNumberOfVisitsPerVisitDuration( $idSite, $period, $date )
	{
		Piwik::checkUserHasViewAccess( $idSite );
		$archive = Piwik_Archive::build($idSite, $period, $date );
		$dataTable = $archive->getDataTable('VisitorInterest_timeGap');
		
		$dataTable->queueFilter('Piwik_DataTable_Filter_ReplaceColumnNames');
		$dataTable->queueFilter('Piwik_DataTable_Filter_Sort', array('label', 'asc', true));
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackReplace', array('label', 'Piwik_getDurationLabel'));
	
		return $dataTable;
	}

	public function getNumberOfVisitsPerPage( $idSite, $period, $date )
	{
		Piwik::checkUserHasViewAccess( $idSite );
		$archive = Piwik_Archive::build($idSite, $period, $date );
		$dataTable = $archive->getDataTable('VisitorInterest_pageGap');
		$dataTable->queueFilter('Piwik_DataTable_Filter_ReplaceColumnNames');
		$dataTable->queueFilter('Piwik_DataTable_Filter_Sort', array('label', 'asc', true));
		$dataTable->queueFilter('Piwik_DataTable_Filter_ColumnCallbackReplace', array('label', 'Piwik_getPageGapLabel'));
	
		return $dataTable;
	}
}

function Piwik_getDurationLabel($label)
{ 
	if(($pos = strpos($label,'-')) !== false)
	{
		$min = substr($label, 0, $pos);
		$max = substr($label, $pos+1);
		
		if($min == 0 || $min == 30)
		{
			$XYSeconds = Piwik_Translate('VisitorInterest_BetweenXYSeconds');
			return sprintf($XYSeconds, $min, $max);
		}
		else
		{
			$min = $min / 60;
			$max = $max / 60;
			$XYMin = Piwik_Translate('VisitorInterest_BetweenXYMinutes');
			return sprintf($XYMin, $min, $max);
		}
	}
	else
	{
		$time = intval($label) / 60;
		$plusXMin = Piwik_Translate('VisitorInterest_PlusXMin');
		return sprintf($plusXMin, urlencode('+').$time);
	}
}

function Piwik_getPageGapLabel($label)
{
	$return = false;
	if(($pos = strpos($label,'-')) !== false)
	{
		$min = substr($label, 0, $pos);
		$max = substr($label, $pos+1);
		
		if($min == $max)
		{
			$return = $min;
		}
	}
	if(!$return)
	{
		$return = $label;
	}
	
	if($return == 1)
	{
		return Piwik_Translate('VisitorInterest_OnePage');
	}
	return sprintf(Piwik_Translate('VisitorInterest_NPages'), $return);
}
