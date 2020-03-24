#!/usr/bin/python

import unittest, shutil, argparse, os
import scans

class TestCommandHelp(unittest.TestCase):
    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in scans.__commands__:
            parser = parser_fun(argparse.ArgumentParser())
            helpstring = parser.format_help()

class TestScans(unittest.TestCase):
    # to be executed prior to running tests
    def setUp(self):
        self.cms_selection = scans.Selection()
        pass

    # to be executed after running tests, regardless of pass/fail
    # only called if setUp succeeds
    def tearDown(self):
        pass

    def test_selection_return_code(self):
        self.assertEqual( self.cms_selection.main(), 0 )
