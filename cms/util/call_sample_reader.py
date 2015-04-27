#!/usr/bin/python

__author__ = "tomkinsc@broadinstitute.org"
__version__ = "PLACEHOLDER"
__date__ = "PLACEHOLDER"

import csv
import os, re
from collections import defaultdict
from urllib2 import urlopen

class CallSampleReader(object):
    '''
        CallSampleReader reads in a tab-delimited file of the sort used to represent 1000 Genomes 
        metadata relating sample IDs (of the sort present in the 1000G VCF files) to
        specific populations, super populations, and gender. It can then filter sample IDs.

        inFile can either be a local file path, or a url

        The input file is expected to have four columns, in any order, with a header described as follows:
        sample  pop super_pop   gender      
        
        An example row:
        HG00096 GBR EUR male

        A fill example input file can be found on the 1000G server:
        ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
    '''
    def __init__(self, inFile):
        self.inFile = inFile

        self.sample_names     = []
        self.sample_membership = defaultdict(dict)
        self.pop_membership    = defaultdict(list)

        self.readFileAndParseWithFunctions(self.store_sample_membership, self.store_population_membership, self.store_population_names)

    @staticmethod
    def read_sample_membership(inFile):
        '''
            Returns a dict where each entry is sample:{super_pop:v, gender:v, pop:v} 
        '''
        sample_membership = defaultdict(dict)
        isLocal = True
        try:
            f = urlopen(inFile)
            isLocal = False
        except ValueError:  # invalid URL
            f = open(inFile, mode='r')

        reader = csv.DictReader(f, delimiter="\t")
        for row in  reader:
            sample_membership[row["sample"]] = {k.replace(" ", "_"):v for k,v in row.iteritems() if k not in ("sample", "")}
        if isLocal:
            f.close()
        return sample_membership

    @staticmethod
    def read_population_membership(inFile):
        '''
            Returns a dict where each entry is pop:[sample, sample, ...]
            The key "ALL" contains all sample names.
        '''
        pop_membership = defaultdict(list)
        isLocal = True
        try:
            f = urlopen(inFile)
            isLocal = False
        except ValueError:  # invalid URL
            f = open(inFile, mode='r')

        reader = csv.DictReader(f, delimiter="\t")
        for row in  reader:
            for field in reader.fieldnames:
                if field not in ("sample",""):
                    pop_membership[row[field]].append(row["sample"])
            pop_membership["ALL"].append(row["sample"])
        if isLocal:
            f.close()
        return pop_membership

    @staticmethod
    def read_sample_names(inFile):
        '''
            Returns a list of sample names
        '''
        sample_names = []
        isLocal = True
        try:
            f = urlopen(inFile)
            isLocal = False
        except ValueError:  # invalid URL
            f = open(inFile, mode='r')

        reader = csv.DictReader(f, delimiter="\t")

        for row in  reader:
            sample_names.append(row["sample"])
        if isLocal:
            f.close()
        return sample_names

    def store_sample_membership(self, row, rowFieldNames):
        '''
            Add the metadata keys for a given row to the sample key of the sample_membership dict
        '''
        self.sample_membership[row["sample"]] = {k.replace(" ", "_"):v for k,v in row.iteritems() if k not in ("sample", "")}


    def store_population_membership(self, row, rowFieldNames):
        '''
            Add the sample to the corresponding pop and super_pop lists
        '''
        for field in rowFieldNames:
            if field not in ("sample",""):
                self.pop_membership[row[field]].append(row["sample"])
        self.pop_membership["ALL"].append(row["sample"])


    def store_population_names(self, row, rowFieldNames):
        '''
            Add the sample to the list of all sample names
        '''
        self.sample_names.append(row["sample"])

    def readFileAndParseWithFunctions(self, *args):
        isLocal = True
        try:
            f = urlopen(self.inFile)
            isLocal = False
        except ValueError:  # invalid URL
            f = open(self.inFile, mode='r')

        reader = csv.DictReader(f, delimiter="\t")
        for row in  reader:
            for func in args:
                func(row, reader.fieldnames)        

        if isLocal:
            f.close()

    def filterSamples(self, **kwargs):
        '''
            Filter the samples to include those matching the specified keywords
            Match patterns can be passed in as either single strings or lists of strings

            The generator returned includes samples matching ANY keyword specified

            For an instance of CallSampleReader(), csr, this function may be used as follows:
            csr.filterSamples(pop="GBR")
            csr.filterSamples(pop=["GBR","FIN"])
            csr.filterSamples(pop=["GBR","FIN"],super_pop="EUR")

            For exclusion (for example all of the super_pop "EUR" except "FIN" and "GBR") use set()
            all_eur_but_fin = list(set(csr.filterSamples(super_pop="EUR")) - set(csr.filterSamples(pop=["FIN"])))
        '''
        for key, value in kwargs.items():
            patterns = []
            # handle passing in either a single string or a list of strings
            if isinstance(value, list):
                patterns = value
            else:
                patterns.append(value)

            return (k for k, v in self.sample_membership.items() if v[key] in patterns)

if __name__ == "__main__":
    rootDir = os.path.dirname(os.path.realpath(__file__))
    filePath = rootDir + "/integrated_call_samples_v3.20130502.ALL.panel"
    fileUrl = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"


    csr = CallSampleReader(filePath)
    #csr = CallSampleReader(fileUrl)

    #print(CallSampleReader.read_population_membership(fileUrl))
    # print("EUR count: ", len( list(csr.filterSamples(super_pop="EUR")) ))
    # print("FIN count: ", len( list(csr.filterSamples(pop=["FIN"])) ))
    # print("FIN count: ", len( list(csr.filterSamples(pop=["GBR"])) ))
    # print("Difference count (EUR-(FIN+GBR)): ", len(list(set(csr.filterSamples(super_pop="EUR")) - set(csr.filterSamples(pop=["FIN","GBR"])))))


