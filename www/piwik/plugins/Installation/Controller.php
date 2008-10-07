<?php
/**
 * Piwik - Open source web analytics
 * 
 * @link http://piwik.org
 * @license http://www.gnu.org/licenses/gpl-3.0.html Gpl v3 or later
 * @version $Id: Controller.php 565 2008-07-21 00:31:13Z matt $
 * 
 * @package Piwik_Installation
 */

require_once "View.php";
require_once "Installation/View.php";

/**
 * 
 * @package Piwik_Installation
 */
class Piwik_Installation_Controller extends Piwik_Controller
{
	// public so plugins can modify it
	public $steps = array(
			'welcome',
			'systemCheck',
			'databaseSetup',
			'tablesCreation',
			'generalSetup',
			'firstWebsiteSetup',
			'displayJavascriptCode',
			'finished'
		);
		
	protected $pathView = 'Installation/templates/';
	
	public function __construct()
	{
		if(!isset($_SESSION['currentStepDone'])) 
		{
			$_SESSION['currentStepDone'] = '';
		}
		
		Piwik_PostEvent('InstallationController.construct', $this);
	}
	
	public function getInstallationSteps()
	{
		return $this->steps;
	}
	
	function getDefaultAction()
	{
		return $this->steps[0];
	}
	
	function welcome()
	{
		$view = new Piwik_Install_View(
						$this->pathView . 'welcome.tpl', 
						$this->getInstallationSteps(),
						__FUNCTION__
					);
		$this->skipThisStep( __FUNCTION__ );
		$view->showNextStep = true;
		
		$_SESSION['currentStepDone'] = __FUNCTION__;		
		echo $view->render();
	}
	
	function systemCheck()
	{
		$view = new Piwik_Install_View(
						$this->pathView . 'systemCheck.tpl', 
						$this->getInstallationSteps(),
						__FUNCTION__
					);
		$this->checkPreviousStepIsValid( __FUNCTION__ );
		$this->skipThisStep( __FUNCTION__ );
		
		$view->infos = $this->getSystemInformation();
		$view->problemWithSomeDirectories = (false !== array_search(false, $view->infos['directories']));
		
		$view->showNextStep = !$view->problemWithSomeDirectories 
							&& $view->infos['phpVersion_ok']
							&& $view->infos['pdo_ok']
							&& $view->infos['pdo_mysql_ok']

						;
		$_SESSION['currentStepDone'] = __FUNCTION__;		

		echo $view->render();
	}
	
	
	function databaseSetup()
	{
		// case the user hits the back button
		$_SESSION['skipThisStep']['firstWebsiteSetup'] = false;
		$_SESSION['skipThisStep']['displayJavascriptCode'] = false;
		
		$view = new Piwik_Install_View(
						$this->pathView . 'databaseSetup.tpl', 
						$this->getInstallationSteps(),
						__FUNCTION__
					);
		$this->checkPreviousStepIsValid( __FUNCTION__ );
		$this->skipThisStep( __FUNCTION__ );
					
		$view->showNextStep = false;
		require_once "FormDatabaseSetup.php";
		$form = new Piwik_Installation_FormDatabaseSetup;
		
		if($form->validate())
		{
			$dbInfos = array(
				'host' 			=> $form->getSubmitValue('host'),
				'username' 		=> $form->getSubmitValue('username'),
				'password' 		=> $form->getSubmitValue('password'),
				'dbname' 		=> $form->getSubmitValue('dbname'),
				'tables_prefix' => $form->getSubmitValue('tables_prefix'),
				'adapter' 		=> Zend_Registry::get('config')->database->adapter,
				'port'			=> Zend_Registry::get('config')->database->port,
			);
			
			try{ 
				$dbInfos['password'] = '"'.htmlspecialchars($form->getSubmitValue('password')).'"';
				
				if(($portIndex = strpos($dbInfos['host'],':')) !== false)
				{
					$dbInfos['port'] = substr($dbInfos['host'], $portIndex + 1 );
					$dbInfos['host'] = substr($dbInfos['host'], 0, $portIndex);
				}
				Piwik::createDatabaseObject($dbInfos);
				
				$_SESSION['db_infos'] = $dbInfos;
				
				$this->redirectToNextStep( __FUNCTION__ );
			} catch(Exception $e) {
				$view->errorMessage = $e->getMessage();
			}
		}
		$view->addForm($form);
		
		$view->infos = $this->getSystemInformation();
			
		echo $view->render();
	}
	
	protected function skipThisStep( $step )
	{
		if(isset($_SESSION['skipThisStep'][$step])
			&& $_SESSION['skipThisStep'][$step])
		{
			$this->redirectToNextStep($step);
		}
	}
	
	function tablesCreation()
	{
		$view = new Piwik_Install_View(
						$this->pathView . 'tablesCreation.tpl', 
						$this->getInstallationSteps(),
						__FUNCTION__
					);
		$this->checkPreviousStepIsValid( __FUNCTION__ );
		$this->skipThisStep( __FUNCTION__ );
		$this->createDbFromSessionInformation();
		
		if(Piwik_Common::getRequestVar('deleteTables', 0, 'int') == 1)
		{
			Piwik::dropTables();
			$view->existingTablesDeleted = true;
			
			// when the user decides to drop the tables then we dont skip the next steps anymore
			$_SESSION['skipThisStep']['firstWebsiteSetup'] = false;
			$_SESSION['skipThisStep']['displayJavascriptCode'] = false;
		}
		
		$tablesInstalled = Piwik::getTablesInstalled();
		$tablesToInstall = Piwik::getTablesNames();
		
		if(count($tablesInstalled) > 0)
		{
			$view->someTablesInstalled = true;
			$view->tablesInstalled = implode(", ", $tablesInstalled);
			
			// when the user reuses the same tables we skip the website creation step
			$_SESSION['skipThisStep']['firstWebsiteSetup'] = true;
			$_SESSION['skipThisStep']['displayJavascriptCode'] = true;
		}
		else
		{
			Piwik::createTables();
			Piwik::createAnonymousUser();
			Piwik::createTablesIndex();
			
			$view->tablesCreated = true;
			$view->showNextStep = true;
		}
		
		$_SESSION['currentStepDone'] = __FUNCTION__;
		echo $view->render();
	}
	
	function generalSetup()
	{		
		$view = new Piwik_Install_View(
						$this->pathView . 'generalSetup.tpl', 
						$this->getInstallationSteps(),
						__FUNCTION__
					);
		$this->checkPreviousStepIsValid( __FUNCTION__ );
		$this->skipThisStep( __FUNCTION__ );
		
		require_once "FormGeneralSetup.php";
		$form = new Piwik_Installation_FormGeneralSetup;
		
		if($form->validate())
		{			
			$superUserInfos = array(
				'login' 		=> $form->getSubmitValue('login'),
				'password' 		=> md5( $form->getSubmitValue('password') ),
				'email' 		=> $form->getSubmitValue('email'),
			);
			
			$_SESSION['superuser_infos'] = $superUserInfos;
			$this->redirectToNextStep( __FUNCTION__ );
		}
		$view->addForm($form);
			
		echo $view->render();
	}
	
	public function firstWebsiteSetup()
	{
				
		$view = new Piwik_Install_View(
						$this->pathView . 'firstWebsiteSetup.tpl', 
						$this->getInstallationSteps(),
						__FUNCTION__
					);
		$this->checkPreviousStepIsValid( __FUNCTION__ );
		$this->skipThisStep( __FUNCTION__ );
		
		require_once "FormFirstWebsiteSetup.php";
		$form = new Piwik_Installation_FormFirstWebsiteSetup;
		
		if( !isset($_SESSION['generalSetupSuccessMessage']))
		{
			$view->displayGeneralSetupSuccess = true;
			$_SESSION['generalSetupSuccessMessage'] = true;
		}
		
		if($form->validate())
		{
			// we setup the superuser login & password in the config that will be checked by the
			// API authentication process
			Zend_Registry::get('config')->superuser = $_SESSION['superuser_infos'];
			
			$name = urlencode($form->getSubmitValue('siteName'));
			$url = urlencode($form->getSubmitValue('url'));
			
			$this->initObjectsToCallAPI();
						
			require_once "API/Request.php";
			$request = new Piwik_API_Request("
							method=SitesManager.addSite
							&siteName=$name
							&urls=$url
							&format=original
						");
						
			try {
				$result = $request->process();
				$_SESSION['site_idSite'] = $result;
				$_SESSION['site_name'] = $name;
				$_SESSION['site_url'] = $url;
				
				$this->redirectToNextStep( __FUNCTION__ );
			} catch(Exception $e) {
				$view->errorMessage = $e->getMessage();
			}

		}
		$view->addForm($form);
		
		echo $view->render();
	}
	public function displayJavascriptCode()
	{
		$view = new Piwik_Install_View(
						$this->pathView . 'displayJavascriptCode.tpl', 
						$this->getInstallationSteps(),
						__FUNCTION__
					);
		$this->checkPreviousStepIsValid( __FUNCTION__ );
		$this->skipThisStep( __FUNCTION__ );
		
		if( !isset($_SESSION['firstWebsiteSetupSuccessMessage']))
		{
			$view->displayfirstWebsiteSetupSuccess = true;
			$_SESSION['firstWebsiteSetupSuccessMessage'] = true;
		}
		
		
		$view->websiteName = urldecode($_SESSION['site_name']);
		
		$jsTag = Piwik::getJavascriptCode($_SESSION['site_idSite'], Piwik_Url::getCurrentUrlWithoutFileName());
		
		$view->javascriptTag = $jsTag;
		$view->showNextStep = true;
		
		$_SESSION['currentStepDone'] = __FUNCTION__;
		echo $view->render();
	}
	
	public function finished()
	{
		$view = new Piwik_Install_View(
						$this->pathView . 'finished.tpl', 
						$this->getInstallationSteps(),
						__FUNCTION__
					);
		$this->writeConfigFileFromSession();
		$this->checkPreviousStepIsValid( __FUNCTION__ );
		$this->skipThisStep( __FUNCTION__ );
		
		
		$_SESSION['currentStepDone'] = __FUNCTION__;		
		$view->showNextStep = false;
		
	    setcookie(session_name(), session_id(), 1, '/');
		@session_destroy();	
		echo $view->render();
		
		// cron tab help
		// javascript reminder
		// giving good names to pages
	}
	
	protected function initObjectsToCallAPI()
	{
		// connect to the database using the DB infos currently in the session
		$this->createDbFromSessionInformation();

		// create the fake access to grant super user privilege
		Zend_Registry::set('access', new Piwik_FakeAccess_SetSuperUser);
		
		// we need to create the logs otherwise the API request throws an exception
		Piwik::createLogObject();
	}
	
	protected function writeConfigFileFromSession()
	{
		$configFile = "; <?php exit; ?> DO NOT REMOVE THIS LINE\n";
		$configFile .= "; file automatically generated during the piwik installation process\n";
		
		// super user information
		$configFile .= "[superuser]\n";
		foreach( $_SESSION['superuser_infos'] as $key => $value)
		{
			$configFile .= "$key = $value\n";
		}
		$configFile .= "\n";
		
		// database information
		$configFile .= "[database]\n";
		foreach($_SESSION['db_infos'] as  $key => $value)
		{
			$configFile .= "$key = $value\n";
		}
		
		file_put_contents(Piwik_Config::getDefaultUserConfigPath(), $configFile);
	}
	/**
	 * The previous step is valid if it is either 
	 * - any step before (OK to go back)
	 * - the current step (case when validating a form)
	 */
	function checkPreviousStepIsValid( $currentStep )
	{
		// the currentStep
		$currentStepId = array_search($currentStep, $this->steps);
		
		// the step before
		$previousStepId = array_search($_SESSION['currentStepDone'], $this->steps);
		
		// not OK if currentStepId > previous+1
		if( $currentStepId > $previousStepId + 1 )
		{
			$message = "Error: it seems you try to skip a step of the Installation process, #
						or your cookies are disabled. 
						<br /><b>Make sure your cookies are enabled</b> and go back 
						<a href='".Piwik_Url::getCurrentUrlWithoutFileName()."'>
						to the first page of the installation</a>.";
			Piwik::exitWithErrorMessage( $message );
		}		
	}

	protected function redirectToNextStep($currentStep)
	{
		$_SESSION['currentStepDone'] = $currentStep;
		$nextStep = $this->steps[1 + array_search($currentStep, $this->steps)];
		Piwik::redirectToModule('Installation' , $nextStep);
	}
	
	protected function createDbFromSessionInformation()
	{
		$dbInfos = $_SESSION['db_infos'];
		
		Zend_Registry::get('config')->database = $dbInfos;
		Piwik::createDatabaseObject($dbInfos);
	}
	
	protected function getSystemInformation()
	{
		$minimumPhpVersion = Zend_Registry::get('config')->General->minimum_php_version;
		$minimumMemoryLimit = Zend_Registry::get('config')->General->minimum_memory_limit;
		
		$infos = array();
	
		$infos['directories'] = Piwik::checkDirectoriesWritable();
		$infos['phpVersion_minimum'] = $minimumPhpVersion;
		$infos['phpVersion'] = phpversion();
		$infos['phpVersion_ok'] = version_compare( $minimumPhpVersion, $infos['phpVersion']) === -1;
		
		$extensions = @get_loaded_extensions();
		
		$infos['pdo_ok'] = false;
		if (in_array('PDO', $extensions))  
		{
		    $infos['pdo_ok'] = true;
		}
				
		$infos['pdo_mysql_ok'] = false;
		if (in_array('pdo_mysql', $extensions))  
		{
		    $infos['pdo_mysql_ok'] = true;
		}
		
		$infos['gd_ok'] = false;
		if (in_array('gd', $extensions)) 
		{
		    $gdInfo = gd_info();
			$infos['gd_version'] = $gdInfo['GD Version'];
		    ereg ("([0-9]{1})", $gdInfo['GD Version'], $gdVersion);
		    if($gdVersion[0] >= 2) 
		    {
				$infos['gd_ok'] = true;
		    }
		}
			
		$infos['serverVersion'] = addslashes($_SERVER['SERVER_SOFTWARE']);
		$infos['serverOs'] = @php_uname();
		$infos['serverTime'] = date('H:i:s');

		$infos['setTimeLimit_ok'] = false;
		if(function_exists( 'set_time_limit'))
		{
			$infos['setTimeLimit_ok'] = true;
		}

		$infos['mail_ok'] = false;
		if(function_exists('mail'))
		{
			$infos['mail_ok'] = true;
		}
		
		$infos['registerGlobals_ok'] = ini_get('register_globals') == 0;
		$infos['memoryMinimum'] = $minimumMemoryLimit;
		
		$infos['memory_ok'] = true;
		// on windows the ini_get is not working?
		$infos['memoryCurrent'] = '?M';

		$raised = Piwik::raiseMemoryLimitIfNecessary();
		if(	$memoryValue = Piwik::getMemoryLimitValue() )
		{
			$infos['memoryCurrent'] = $memoryValue."M";
			$infos['memory_ok'] = $memoryValue >= $minimumMemoryLimit;
		}
		
		return $infos;
	}
}


/**
 * 
 * @package Piwik_Installation
 */
class Piwik_FakeAccess_SetSuperUser {
	function checkUserIsSuperUser()
	{
		return true;
	}
	function loadAccess() {}
}
