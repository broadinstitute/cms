#!/usr/bin/python

import unittest
import main

class TestCMS(unittest.TestCase):
    # to be executed prior to running tests
    def setUp(self):
        self.cms = main.CMS()
        pass

    # to be executed after running tests, regardless of pass/fail
    # only called if setUp succeeds
    def tearDown(self):
        pass

    def test_main_return_code(self):
        self.assertEqual( self.cms.main(), 0 )
