#!/usr/bin/python

import unittest
from selection import Selection

class TestSelection(unittest.TestCase):
    # to be executed prior to running tests
    def setUp(self):
        self.cms_selection = Selection()
        pass

    # to be executed after running tests, regardless of pass/fail
    # only called if setUp succeeds
    def tearDown(self):
        pass

    def test_selection_return_code(self):
        self.assertEqual( self.cms_selection.main(), 0 )
