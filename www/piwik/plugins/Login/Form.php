<?php
/**
 * Piwik - Open source web analytics
 * 
 * @link http://piwik.org
 * @license http://www.gnu.org/licenses/gpl-3.0.html Gpl v3 or later
 * @version $Id: Form.php 550 2008-07-05 01:03:04Z matt $
 * 
 * @package Piwik_Login
 */

require_once "modules/Form.php";

/**
 * 
 * @package Piwik_Login
 */
class Piwik_Login_Form extends Piwik_Form
{
	function __construct()
	{
		parent::__construct();
		// reset 
		$this->updateAttributes('id="loginform" name="loginform"');
	}
	
	function init()
	{
		// if form_url is not defined go to current url
		$urlToGoAfter = Piwik_Common::getRequestVar('form_url', Piwik_Url::getCurrentUrl(), 'string');
		
		// if the current url to redirect contains module=login we insteaed redirect to the referer url
		if(stripos($urlToGoAfter,'module=login') !== false)
		{
			$urlToGoAfter = Piwik_Url::getReferer();
		}
		
		$formElements = array(
			array('text', 'form_login'),
			array('password', 'form_password'),
			array('hidden', 'form_url', $urlToGoAfter),
		);
		$this->addElements( $formElements );
		
		$formRules = array(
			array('form_login', sprintf(Piwik_Translate('General_Required'), Piwik_Translate('Login_Login')), 'required'),
			array('form_password', sprintf(Piwik_Translate('General_Required'), Piwik_Translate('Login_Password')), 'required'),
		);
		$this->addRules( $formRules );	
		$this->addElement('submit', 'submit');
	}
}

