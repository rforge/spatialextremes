<?php
/**
 * Piwik - Open source web analytics
 * 
 * @link http://piwik.org
 * @license http://www.gnu.org/licenses/gpl-3.0.html Gpl v3 or later
 * @version $Id: Referers.php 561 2008-07-21 00:00:35Z matt $
 * 
 * @package Piwik_Referers
 */
	
/**
 * 
 * @package Piwik_Referers
 */
class Piwik_Referers extends Piwik_Plugin
{	
	public function setCategoryDelimiter($delimiter)
	{
		self::$actionCategoryDelimiter = $delimiter;
	}
	

	public function getInformation()
	{
		$info = array(
			'name' => 'Referers',
			'description' => 'Computes all reports about the referers',
			'author' => 'Piwik',
			'homepage' => 'http://piwik.org/',
			'version' => '0.1',
		);
		
		return $info;
	}
	
	function postLoad()
	{
		Piwik_AddWidget( 'Referers', 'getKeywords', Piwik_Translate('Referers_WidgetKeywords'));
		Piwik_AddWidget( 'Referers', 'getPartners', Piwik_Translate('Referers_WidgetPartners'));
		Piwik_AddWidget( 'Referers', 'getCampaigns', Piwik_Translate('Referers_WidgetCampaigns'));
		Piwik_AddWidget( 'Referers', 'getWebsites', Piwik_Translate('Referers_WidgetExternalWebsites'));
		Piwik_AddWidget( 'Referers', 'getSearchEngines', Piwik_Translate('Referers_WidgetSearchEngines'));
		Piwik_AddWidget( 'Referers', 'getRefererType', Piwik_Translate('Referers_WidgetOverview'));

		Piwik_AddMenu('Referers', Piwik_Translate('Referers_SubmenuEvolution'), array('module' => 'Referers'));
		Piwik_AddMenu('Referers', Piwik_Translate('Referers_SubmenuSearchEngines'), array('module' => 'Referers', 'action' => 'getSearchEnginesAndKeywords'));
		Piwik_AddMenu('Referers', Piwik_Translate('Referers_SubmenuWebsites'), array('module' => 'Referers', 'action' => 'getWebsites'));
		Piwik_AddMenu('Referers', Piwik_Translate('Referers_SubmenuCampaigns'), array('module' => 'Referers', 'action' => 'getCampaigns'));
		Piwik_AddMenu('Referers', Piwik_Translate('Referers_SubmenuPartners'), array('module' => 'Referers', 'action' => 'getPartners'));	
	}
		
	function getListHooksRegistered()
	{
		$hooks = array(
			'ArchiveProcessing_Day.compute' => 'archiveDay',
			'ArchiveProcessing_Period.compute' => 'archiveMonth',
		);
		return $hooks;
	}
	
	
	function archiveMonth( $notification )
	{
		$archiveProcessing = $notification->getNotificationObject();
		
		$dataTableToSum = array( 
				'Referers_type',
				'Referers_keywordBySearchEngine',
				'Referers_searchEngineByKeyword',
				'Referers_keywordByCampaign',
				'Referers_urlByWebsite',
				'Referers_urlByPartner',
		);
		
		$maximumRowsInDataTableLevelZero = 500;
		$maximumRowsInSubDataTable = 50;
		$nameToCount = $archiveProcessing->archiveDataTable($dataTableToSum, $maximumRowsInDataTableLevelZero, $maximumRowsInSubDataTable);
		
		$mappingFromArchiveName = array(
			'Referers_distinctSearchEngines' => 
						array( 	'typeCountToUse' => 'level0',
								'nameTableToUse' => 'Referers_keywordBySearchEngine',
							),
			'Referers_distinctKeywords' => 
						array( 	'typeCountToUse' => 'level0',
								'nameTableToUse' => 'Referers_searchEngineByKeyword',
							),
			'Referers_distinctCampaigns' => 
						array( 	'typeCountToUse' => 'level0',
								'nameTableToUse' => 'Referers_keywordByCampaign',
							),
			'Referers_distinctWebsites' => 
						array( 	'typeCountToUse' => 'level0',
								'nameTableToUse' => 'Referers_urlByWebsite',
							),
			'Referers_distinctWebsitesUrls' => 
						array( 	'typeCountToUse' => 'recursive',
								'nameTableToUse' => 'Referers_urlByWebsite',
							),
			'Referers_distinctPartners' => 
						array( 	'typeCountToUse' => 'level0',
								'nameTableToUse' => 'Referers_urlByPartner',
							),
			'Referers_distinctPartnersUrls' => 
						array( 	'typeCountToUse' => 'recursive',
								'nameTableToUse' => 'Referers_urlByPartner',
							),
		);

		foreach($mappingFromArchiveName as $name => $infoMapping)
		{
			$typeCountToUse = $infoMapping['typeCountToUse'];
			$nameTableToUse = $infoMapping['nameTableToUse'];
			
			if($typeCountToUse == 'recursive')
			{
				
				$countValue = $nameToCount[$nameTableToUse]['recursive']
								- $nameToCount[$nameTableToUse]['level0'];
			}
			else
			{
				$countValue = $nameToCount[$nameTableToUse]['level0'];
			}
			
			$record = new Piwik_ArchiveProcessing_Record_Numeric(
													$name, 
													$countValue
												);
		}
	}

	public function archiveDay( $notification )
	{
		$archiveProcessing = $notification->getNotificationObject();
		$query = "SELECT 	referer_type, 
							referer_name, 
							referer_keyword,
							referer_url,
							count(distinct visitor_idcookie) as nb_uniq_visitors,
							count(*) as nb_visits,
							sum(visit_total_actions) as nb_actions,
							max(visit_total_actions) as max_actions, 
							sum(visit_total_time) as sum_visit_length,							
							sum(case visit_total_actions when 1 then 1 else 0 end) as bounce_count
				 	FROM ".$archiveProcessing->logTable."
				 	WHERE visit_server_date = ?
				 		AND idsite = ?
				 	GROUP BY referer_type, referer_name, referer_keyword";
		$query = $archiveProcessing->db->query($query, array( $archiveProcessing->strDateStart, $archiveProcessing->idsite ));
		
		$timer = new Piwik_Timer;
		
		
		$interestBySearchEngine =
			$interestByKeyword =
			$keywordBySearchEngine =
			$searchEngineByKeyword =
			$interestByWebsite[Piwik_Common::REFERER_TYPE_WEBSITE] =
			$interestByWebsite[Piwik_Common::REFERER_TYPE_PARTNER] =
			$urlByWebsite[Piwik_Common::REFERER_TYPE_WEBSITE] =
			$urlByWebsite[Piwik_Common::REFERER_TYPE_PARTNER] =
			$interestByNewsletter =
			$keywordByCampaign =
			$interestByCampaign =
			$interestByType = 
			$distinctUrls[Piwik_Common::REFERER_TYPE_WEBSITE] =
			$distinctUrls[Piwik_Common::REFERER_TYPE_PARTNER] = array();
		
		while($rowBefore = $query->fetch() )
		{
			$row = array(
				Piwik_Archive::INDEX_NB_UNIQ_VISITORS 	=> $rowBefore['nb_uniq_visitors'], 
				Piwik_Archive::INDEX_NB_VISITS 			=> $rowBefore['nb_visits'], 
				Piwik_Archive::INDEX_NB_ACTIONS 		=> $rowBefore['nb_actions'], 
				Piwik_Archive::INDEX_MAX_ACTIONS 		=> $rowBefore['max_actions'], 
				Piwik_Archive::INDEX_SUM_VISIT_LENGTH 	=> $rowBefore['sum_visit_length'], 
				Piwik_Archive::INDEX_BOUNCE_COUNT 		=> $rowBefore['bounce_count'],
				'referer_type' 							=> $rowBefore['referer_type'],
				'referer_name' 							=> $rowBefore['referer_name'], 
				'referer_keyword'						=> $rowBefore['referer_keyword'],
				'referer_url'							=> $rowBefore['referer_url'],
			);
			
			switch($row['referer_type'])
			{
				case Piwik_Common::REFERER_TYPE_SEARCH_ENGINE:
				
					if(!isset($interestBySearchEngine[$row['referer_name']])) $interestBySearchEngine[$row['referer_name']]= $archiveProcessing->getNewInterestRow();
					if(!isset($interestByKeyword[$row['referer_keyword']])) $interestByKeyword[$row['referer_keyword']]= $archiveProcessing->getNewInterestRow();
					if(!isset($keywordBySearchEngine[$row['referer_name']][$row['referer_keyword']])) $keywordBySearchEngine[$row['referer_name']][$row['referer_keyword']]= $archiveProcessing->getNewInterestRow();
					if(!isset($searchEngineByKeyword[$row['referer_keyword']][$row['referer_name']])) $searchEngineByKeyword[$row['referer_keyword']][$row['referer_name']]= $archiveProcessing->getNewInterestRow();
				
					$archiveProcessing->updateInterestStats( $row, $interestBySearchEngine[$row['referer_name']]);
					$archiveProcessing->updateInterestStats( $row, $interestByKeyword[$row['referer_keyword']]);
					$archiveProcessing->updateInterestStats( $row, $keywordBySearchEngine[$row['referer_name']][$row['referer_keyword']]);
					$archiveProcessing->updateInterestStats( $row, $searchEngineByKeyword[$row['referer_keyword']][$row['referer_name']]);
				break;
				
				case Piwik_Common::REFERER_TYPE_WEBSITE:
				case Piwik_Common::REFERER_TYPE_PARTNER:
					
					if(!isset($interestByWebsite[$row['referer_type']][$row['referer_name']])) $interestByWebsite[$row['referer_type']][$row['referer_name']]= $archiveProcessing->getNewInterestRow();
					$archiveProcessing->updateInterestStats( $row, $interestByWebsite[$row['referer_type']][$row['referer_name']]);
					
					if(!isset($urlByWebsite[$row['referer_type']][$row['referer_name']][$row['referer_url']])) $urlByWebsite[$row['referer_type']][$row['referer_name']][$row['referer_url']]= $archiveProcessing->getNewInterestRow();
					$archiveProcessing->updateInterestStats( $row, $urlByWebsite[$row['referer_type']][$row['referer_name']][$row['referer_url']]);
				
					if(!isset($distinctUrls[$row['referer_type']][$row['referer_url']]))
					{
						$distinctUrls[$row['referer_type']][$row['referer_url']] = true;
					}
					
				break;
				
				case Piwik_Common::REFERER_TYPE_NEWSLETTER:
					if(!isset($interestByNewsletter[$row['referer_name']])) $interestByNewsletter[$row['referer_name']]= $archiveProcessing->getNewInterestRow();
					$archiveProcessing->updateInterestStats( $row, $interestByNewsletter[$row['referer_name']]);
					
				break;
				
				case Piwik_Common::REFERER_TYPE_CAMPAIGN:
					if(!empty($row['referer_keyword']))
					{
						if(!isset($keywordByCampaign[$row['referer_name']][$row['referer_keyword']])) $keywordByCampaign[$row['referer_name']][$row['referer_keyword']]= $archiveProcessing->getNewInterestRow();
						$archiveProcessing->updateInterestStats( $row, $keywordByCampaign[$row['referer_name']][$row['referer_keyword']]);
					}
					if(!isset($interestByCampaign[$row['referer_name']])) $interestByCampaign[$row['referer_name']]= $archiveProcessing->getNewInterestRow();
					$archiveProcessing->updateInterestStats( $row, $interestByCampaign[$row['referer_name']]);
				break;
			}
			
			if(!isset($interestByType[$row['referer_type']] )) $interestByType[$row['referer_type']] = $archiveProcessing->getNewInterestRow();
			$archiveProcessing->updateInterestStats($row, $interestByType[$row['referer_type']]);
		}
		
		$numberOfDistinctSearchEngines = count($keywordBySearchEngine);
		$numberOfDistinctKeywords = count($searchEngineByKeyword);
		
		$numberOfDistinctCampaigns = count($interestByCampaign);
		$numberOfDistinctWebsites = count($interestByWebsite[Piwik_Common::REFERER_TYPE_WEBSITE]);
		$numberOfDistinctWebsitesUrls = count($distinctUrls[Piwik_Common::REFERER_TYPE_WEBSITE]);
		$numberOfDistinctPartners = count($interestByWebsite[Piwik_Common::REFERER_TYPE_PARTNER]);
		$numberOfDistinctPartnersUrls = count($distinctUrls[Piwik_Common::REFERER_TYPE_PARTNER]);
		
		$numericRecords = array(
			'Referers_distinctSearchEngines'	=> $numberOfDistinctSearchEngines,
			'Referers_distinctKeywords' 		=> $numberOfDistinctKeywords,
			'Referers_distinctCampaigns'		=> $numberOfDistinctCampaigns,
			'Referers_distinctWebsites'			=> $numberOfDistinctWebsites,
			'Referers_distinctWebsitesUrls'		=> $numberOfDistinctWebsitesUrls,
			'Referers_distinctPartners'			=> $numberOfDistinctPartners,
			'Referers_distinctPartnersUrls'		=> $numberOfDistinctPartnersUrls,
		);
		foreach($numericRecords as $name => $value)
		{
			$record = new Piwik_ArchiveProcessing_Record_Numeric($name, $value);
		}
		
		$data = $archiveProcessing->getDataTableSerialized($interestByType);
		$record = new Piwik_ArchiveProcessing_Record_BlobArray('Referers_type', $data);
		
		$maximumRowsInDataTableLevelZero = 500;
		$maximumRowsInSubDataTable = 50;
		
		$data = $archiveProcessing->getDataTablesSerialized($keywordBySearchEngine, $interestBySearchEngine, $maximumRowsInDataTableLevelZero, $maximumRowsInSubDataTable);
		$record = new Piwik_ArchiveProcessing_Record_BlobArray('Referers_keywordBySearchEngine', $data);
		
		$data = $archiveProcessing->getDataTablesSerialized($searchEngineByKeyword, $interestByKeyword, $maximumRowsInDataTableLevelZero, $maximumRowsInSubDataTable);
		$record = new Piwik_ArchiveProcessing_Record_BlobArray('Referers_searchEngineByKeyword', $data);
		
		$data = $archiveProcessing->getDataTablesSerialized($keywordByCampaign, $interestByCampaign, $maximumRowsInDataTableLevelZero, $maximumRowsInSubDataTable);
		$record = new Piwik_ArchiveProcessing_Record_BlobArray('Referers_keywordByCampaign', $data);
		
		$data = $archiveProcessing->getDataTablesSerialized($urlByWebsite[Piwik_Common::REFERER_TYPE_WEBSITE], $interestByWebsite[Piwik_Common::REFERER_TYPE_WEBSITE], $maximumRowsInDataTableLevelZero, $maximumRowsInSubDataTable);
		$record = new Piwik_ArchiveProcessing_Record_BlobArray('Referers_urlByWebsite', $data);
		
		$data = $archiveProcessing->getDataTablesSerialized($urlByWebsite[Piwik_Common::REFERER_TYPE_PARTNER], $interestByWebsite[Piwik_Common::REFERER_TYPE_PARTNER], $maximumRowsInDataTableLevelZero, $maximumRowsInSubDataTable);
		$record = new Piwik_ArchiveProcessing_Record_BlobArray('Referers_urlByPartner', $data);
	}
}

