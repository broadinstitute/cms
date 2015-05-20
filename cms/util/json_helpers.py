#!/usr/bin/python

import os
import json

class JSONHelper(object):

    @staticmethod
    def annotate_json_file(file_path, dict_to_set, create_if_nonexistent=True):
        filepath = os.path.abspath(file_path) 

        fileExists = os.path.isfile(filepath)

        metaDataDict = dict_to_set

        if fileExists:
            with open(filepath, "r") as inFile:                
                metaDataDict = json.load(inFile)
            
                for k,v in dict_to_set.iteritems():
                    metaDataDict[k] = v
                    # if None is specified for the value, remove from the dict
                    if v == None:
                        del metaDataDict[k]

        if fileExists or create_if_nonexistent:

            with open(filepath+".temp", "w") as outFile:
                json.dump(metaDataDict, outFile)
        
            # move the newly-annotated file to replace the old one
            os.system('mv {} {}'.format(filepath+".temp", filepath))