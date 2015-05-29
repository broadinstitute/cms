#!/usr/bin/python

import os
import json

class JSONHelper(object):

    @staticmethod
    def annotate_json_file(file_path, dict_to_set, create_if_nonexistent=True, key_to_act_on=None, append=False):
        filepath = os.path.abspath(file_path) 

        fileExists = os.path.isfile(filepath)

        metaDataDict = dict_to_set

        if fileExists:
            with open(filepath, "r") as inFile:                
                metaDataDict = json.load(inFile)

                if append:
                    metaDataDict.setdefault(key_to_act_on, []).append(dict_to_set)
                else:
            
                    for k,v in dict_to_set.iteritems():
                        if not key_to_act_on:
                            metaDataDict[k] = v
                        else:
                            metaDataDict[key_to_act_on] = v
                        # if None is specified for the value, remove from the dict
                        if v == None:
                            del metaDataDict[k]

        if fileExists or create_if_nonexistent:

            with open(filepath+".temp", "w") as outFile:
                json.dump(metaDataDict, outFile, sort_keys=True, indent=4, separators=(',', ': '))
        
            # move the newly-annotated file to replace the old one
            os.system('mv {} {}'.format(filepath+".temp", filepath))