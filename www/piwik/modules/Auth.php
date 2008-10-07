<?php
/**
 * Piwik - Open source web analytics
 *
 * @link http://piwik.org
 * @license http://www.gnu.org/licenses/gpl-3.0.html Gpl v3 or later
 * @version $Id: Auth.php 558 2008-07-20 23:10:38Z matt $
 *
 * @package Piwik
 */

interface Piwik_Auth {
	const SUCCESS_SUPERUSER_AUTH_CODE = 42;
	
	/**
	 * @return string
	 */
	public function getTokenAuth();
	
	/**
	 * @return Piwik_Auth_Result
	 */
	public function authenticate();
}

/**
 *
 * @package Piwik
 */
class Piwik_Auth_Result extends Zend_Auth_Result
{
	public function __construct($code, $identity, array $messages = array())
	{
		// Piwik_Auth::SUCCESS_SUPERUSER_AUTH_CODE, Zend_Auth_Result::SUCCESS, Zend_Auth_Result::FAILURE  
		$this->_code		= (int)$code;
		$this->_identity	= $identity;
		$this->_messages	= $messages;
	}
}
