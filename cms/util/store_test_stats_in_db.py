#!/usr/bin/python

# built-ins
import os
import itertools
from functools import wraps
import json

# external
import peewee as pw
from playhouse import shortcuts as plsc # from pip module "peewee"

class ScanTable(pw.Model):
    pass
    #class Meta:
    #    database = db

class SuperPopulation(ScanTable):
    name = pw.CharField(unique=True)

    def __str__(self):
        return self.name

class Population(ScanTable):
    super_population = pw.ForeignKeyField(SuperPopulation, null=True)

    name = pw.CharField(unique=True)

    def __str__(self):
        return self.name

class LocusInfo(ScanTable):
    '''
        Allele information for a given locus
    '''
    chrom                  = pw.CharField(max_length=32) 
    variant_id             = pw.CharField(max_length=64) 
    pos_bp                 = pw.IntegerField()
    map_pos_cm             = pw.DoubleField()
    ref_allele             = pw.CharField()
    alt_allele             = pw.CharField()
    ancestral_call         = pw.CharField()

    class Meta:
        indexes = (
            (('chrom','pos_bp'), True),
        )

    def __str__(self):
        return ",".join([str(chrom), str(variant_id), str(pos_bp), 
            str(map_pos_cm), str(ref_allele), str(alt_allele), 
            str(ancestral_call)])

class PopulationSet(ScanTable):
    '''
        Class to bundle many populations, with optional name and description.
    '''
    name              = pw.CharField()
    description       = pw.CharField()
    populations       = plsc.ManyToManyField(Population, 
                                                related_name='population_sets')
    super_populations = plsc.ManyToManyField(SuperPopulation, 
                                          related_name='super_population_sets')

    def __str__(self):
        return ",".join([str(pop) for pop in self.populations])

class AlleleInfo(ScanTable):
    '''
        Population-specific allele info
    '''
    locus          = pw.ForeignKeyField(LocusInfo)
    population_set = pw.ForeignKeyField(PopulationSet)

    allele_freq_in_pop_set = pw.DoubleField()

    class Meta:
        indexes = (
            (('locus','population_set'), True),
        )

    def __str__(self):
        return ",".join([str(self.locus), str(self.population_set), 
            str(self.allele_freq_in_pop_set)])

class EHHInfo(ScanTable):
    '''
        EHH scores for a given core locus within a given population
    '''
    population_set = pw.ForeignKeyField(PopulationSet)
    core_locus     = pw.ForeignKeyField(LocusInfo)

    pos_bp_delta     = pw.IntegerField()
    map_pos_cm_delta = pw.DoubleField()
    EHH_ones         = pw.DoubleField()
    EHH_zeros        = pw.DoubleField()

    class Meta:
        indexes = (
            (('population_set','core_locus'), True),
        )

    def __str__(self):
        return ",".join([str(self.core_locus), str(self.population_set), 
            str(self.allele_freq_in_pop_set)])

class IHSInfo(ScanTable):
    '''
       EHH scores for a given locus within a given population 
    '''
    population_set = pw.ForeignKeyField(PopulationSet)
    locus          = pw.ForeignKeyField(LocusInfo)

    IHH_1              = pw.DoubleField()
    IHH_0              = pw.DoubleField()
    IHS_unstandardized = pw.DoubleField()
    IHS_standardized   = pw.DoubleField()

    class Meta:
        indexes = (
            (('locus','population_set'), True),
        )

    def __str__(self):
        return ",".join([str(population_set), str(locus), str(IHH_1), 
            str(IHH_0), str(IHS_unstandardized), str(IHS_standardized)])

class XPEHHInfo(ScanTable):
    '''
        XPEHH scores for a given locus for population A with respect to B
    '''
    locus            = pw.ForeignKeyField(LocusInfo)
    population_set_a = pw.ForeignKeyField(PopulationSet, 
                                                related_name='population_set_a')
    population_set_b = pw.ForeignKeyField(PopulationSet, 
                                                related_name='population_set_b')
    
    population_set_a_one_freq = pw.DoubleField()
    IHH_a                     = pw.DoubleField()
    population_set_b_one_freq = pw.DoubleField()
    IHH_b                     = pw.DoubleField()
    XPEHH_unstandardized      = pw.DoubleField()
    XPEHH_standardized        = pw.DoubleField()

    class Meta:
        indexes = (
            (('locus','population_set_a','population_set_a'), True),
        )

    def __str__(self):
        return ",".join([str(locus), str(population_set_a), 
            str(population_set_b), str(population_set_a_one_freq), str(IHH_a), 
            str(population_set_b_one_freq), str(IHH_b), 
            str(XPEHH_unstandardized), str(XPEHH_standardized)])

# ===========================
# TODO: maybe don't include foreign key values in the string representations
#       of the various tables?

class DatabaseManager(object):
    def __init__(self, db_path):
        PopulationPopulationSet = PopulationSet.populations.get_through_model()
        SuperPopulationPopulationSet = PopulationSet.super_populations.get_through_model()
        self.tables = [ SuperPopulation, 
                            LocusInfo, 
                            Population, 
                            PopulationSet, 
                            AlleleInfo, 
                            EHHInfo, 
                            IHSInfo, 
                            XPEHHInfo,
                            PopulationPopulationSet,
                            SuperPopulationPopulationSet]

        db_path = os.path.abspath( os.path.expanduser(db_path) )
        self.db = pw.SqliteDatabase(db_path)

        # add the db to the table Meta class, rather than replying on peewee's
        # Using construct
        for model in self.tables:
            model._meta.database = self.db

        self.db.connect()
        self.db.create_tables( self.tables, safe=True)
    
    def __del__(self):
        self.db.close()

    # as an alterative to adding the table context to each table
    # we can use the @use_db decorator within this class
    def use_db(f):
        '''
            This wraps functions that make use of the database so they have
            the appropriate database configuration in context.
        '''
        @wraps(f)
        def wrapped(inst, *args, **kwargs):
            # this is peewee's way of setting the db within scope
            with pw.Using(inst.db, inst.tables):
                return f(inst, *args, **kwargs)
        return wrapped
    
    def add_super_population(self, superpop_name):
        item, created = SuperPopulation.get_or_create(name=superpop_name)
        return item

    def add_population(self, pop_name):
        item, created = Population.get_or_create(name=pop_name)
        return item

    def add_locus_info(self, chrom, variant_id, pos_bp, map_pos_cm, ref_allele, 
                        alt_allele, ancestral_call):
        item, created = LocusInfo.get_or_create( chrom          = chrom, 
                                              variant_id     = variant_id, 
                                              pos_bp         = pos_bp, 
                                              map_pos_cm     = map_pos_cm, 
                                              ref_allele     = ref_allele, 
                                              alt_allele     = alt_allele, 
                                              ancestral_call = ancestral_call)
        return item

    def add_population_set(self, population_list, superpopulation_list, 
        name="", description=""):
        pop_set = PopulationSet.create(name=name, description=description)

        for pop_name in population_list:
            pop, created = Population.get_or_create(name=pop_name)
            pop_set.populations.add(pop)

        for super_pop_name in superpopulation_list:
            super_pop, created = SuperPopulation.get_or_create(
                                                        name=super_pop_name)
            pop_set.super_populations.add(super_pop)

        return pop_set

    def add_allele_info(self, locus, population_set, allele_freq_in_pop_set):
        item, created = AlleleInfo.get_or_create(locus=locus, 
                                                population_set=population_set, 
                                allele_freq_in_pop_set=allele_freq_in_pop_set)
        return item

    def add_ehh_info(self, population_set, core_locus, pos_bp_delta, 
        map_pos_cm_delta, EHH_ones, EHH_zeros):
        item, created = EHHInfo.get_or_create( population_set  = population_set, 
                                            core_locus       = core_locus, 
                                            pos_bp_delta     = pos_bp_delta, 
                                            map_pos_cm_delta = map_pos_cm_delta, 
                                            EHH_ones         = EHH_ones, 
                                            EHH_zeros        = EHH_zeros)
        return item

    def add_ihs_info(self, population_set, locus, IHH_1, IHH_0, 
        IHS_unstandardized, IHS_standardized):
        item, created = IHSInfo.get_or_create( population_set = population_set, 
                                    locus               = locus, 
                                    IHH_1               = IHH_1, 
                                    IHH_0               = IHH_0, 
                                    IHS_unstandardized  = IHS_unstandardized, 
                                    IHS_standardized    = IHS_standardized)
        return item

    def add_xpehh_info(self, locus, population_set_a, population_set_b, 
        population_set_a_one_freq, IHH_a, population_set_b_one_freq, IHH_b, 
        XPEHH_standardized, XPEHH_unstandardized):

        item, created = XPEHHInfo.get_or_create( locus = locus,
                        population_set_a          = population_set_a,
                        population_set_b          = population_set_b,
                        population_set_a_one_freq = population_set_a_one_freq,
                        IHH_a                     = IHH_a,
                        population_set_b_one_freq = population_set_b_one_freq,
                        IHH_b                     = IHH_b,
                        XPEHH_standardized        = XPEHH_standardized,
                        XPEHH_unstandardized      = XPEHH_unstandardized)

        return item

    def get_popset_for_pop_names(self, name_list):
        '''
            Get the population set corresponding to the list of population names 
            specified in name_list.
        '''
        pops = Population.select(Population.id).where(Population.name << name_list)

        pop_ids = []
        for pop in pops:
            pop_ids.append(int(pop.id))

        query = PopulationSet.select()

        popSets = {}
        for row in query:
            tempList = []
            for popset in row.populations:
                tempList.append(popset.id)
            popSets[row.id] = set(tempList)

        matchingPopSets = []
        for k,v in popSets.iteritems():
            if v == set(pop_ids):
                matchingPopSets.append(k)

        matchingPopSet = PopulationSet.get(PopulationSet.id==matchingPopSets[0])

        return matchingPopSet

    def get_locus_info(self, chrom, pos_bp):
        return LocusInfo.select().where((LocusInfo.chrom == chrom) 
                                            & (LocusInfo.pos_bp == LocusInfo))

    def get_allele_info(self, locus, population_set):
        return AlleleInfo.select().where( (AlleleInfo.locus.id == locus.id) 
                        & (AlleleInfo.population_set.id == population_set.id) )

    def get_ehh_info(self, core_locus, population_set):
        return EHHInfo.select().where( (EHHInfo.core_locus.id == core_locus.id) 
                            & (EHHInfo.population_set.id == population_set.id) )

    def get_ihs_info(self, locus, population_set):
        return IHSInfo.select().where( (IHSInfo.locus.id == locus.id) 
                            & (IHSInfo.population_set.id == population_set.id) )        

    def get_xpehh_info(self, locus, population_set_a, population_set_b):
        return XPEHHInfo.select().where( (XPEHHInfo.locus.id == locus.id) & 
            (XPEHHInfo.population_set_a.id == population_set_a.id) & 
            (XPEHHInfo.population_set_b.id == population_set_b.id))

class ReadCalculationData(object):
    def __init__(self, jsonMetadataFile):
        self.jsonMetadataFile = jsonMetadataFile
        self.metadata = self.read_metadata( jsonMetadataFile )

    # TODO: move to utils.json_helper
    @classmethod
    def read_metadata(cls, configfile_path):
        abspath = os.path.abspath(configfile_path) 

        fileExists = os.path.isfile(abspath)

        metaDataDict = {}

        if fileExists:
            with open(abspath, "r") as inFile:                
                metaDataDict = json.load(inFile)

        return metaDataDict

if __name__ == "__main__":
    db = DatabaseManager('~/Desktop/scan_stats.db')
    #db.add_locus_info(1, 123, 100, 1.2, "A", "T", "A")
    #result = db.add_population_set( population_list=["FIN","GBR","YRI"], superpopulation_list=["EUR","AFR"] )

    pop_set = db.get_popset_for_pop_names(["FIN","GBR","YRI"])
    print [i.id for i in pop_set.populations]























