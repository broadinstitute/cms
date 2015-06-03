#!/usr/bin/python

# built-ins
import os
import itertools, logging
from functools import wraps
import json
import csv

# module
import util.file, util.json_helpers

# external
import boltons
import peewee as pw
from playhouse import shortcuts as plsc # from pip module "peewee"

log = logging.getLogger(__name__)

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
    ancestral_call         = pw.CharField(null=True)

    class Meta:
        indexes = (
            (('chrom','pos_bp','alt_allele'), True),
        )

    def __str__(self):
        return ",".join([str(self.chrom), str(self.variant_id), str(self.pos_bp), 
            str(self.map_pos_cm), str(self.ref_allele), str(self.alt_allele), 
            str(self.ancestral_call)])

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

    allele_freq_in_pop_set = pw.DoubleField(null=True)

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
    EHH              = pw.DoubleField()

    class Meta:
        indexes = (
            (('population_set','core_locus','pos_bp_delta'), True),
        )

    def __str__(self):
        return ",".join([str(self.core_locus), str(self.population_set), 
            str(self.EHH_ones), str(self.EHH_zeroes), str(self.EHH)])

class IHSInfo(ScanTable):
    '''
       EHH scores for a given locus within a given population 
    '''
    population_set = pw.ForeignKeyField(PopulationSet)
    locus          = pw.ForeignKeyField(LocusInfo)

    ones_freq          = pw.DoubleField()
    IHH_1              = pw.DoubleField()
    IHH_0              = pw.DoubleField()
    IHS_unstandardized = pw.DoubleField()
    IHS_standardized   = pw.DoubleField(null=True)

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
    XPEHH_standardized        = pw.DoubleField(null=True)

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
        SuperPopulationPopulationSet = \
                            PopulationSet.super_populations.get_through_model()
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
        self.db.create_tables( self.tables, safe=True )
    
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
        item, created = LocusInfo.get_or_create( chrom   = chrom, 
                                          variant_id     = variant_id, 
                                          pos_bp         = int(pos_bp),
                                          map_pos_cm     = float(map_pos_cm), 
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
                          allele_freq_in_pop_set=float(allele_freq_in_pop_set))   
        return item

    def add_ehh_info(self, population_set, core_locus, pos_bp_delta, 
        map_pos_cm_delta, EHH_ones, EHH_zeros, EHH):

        try:
            item, created = EHHInfo.get_or_create( 
                                    population_set   = population_set, 
                                    core_locus       = core_locus, 
                                    pos_bp_delta     = int(pos_bp_delta),
                                    map_pos_cm_delta = float(map_pos_cm_delta),
                                    EHH_ones         = float(EHH_ones),
                                    EHH_zeros        = float(EHH_zeros),
                                    EHH              = float(EHH))
            return item
        except pw.IntegrityError:
            log.warning("Second variant considered " +\
            "in EHH calculation at locus {} + {}bp; skipping storage.".format(
                                str(core_locus.variant_id), str(pos_bp_delta)))
            

    def add_ihs_info(self, population_set, locus, ones_freq, IHH_1, IHH_0, 
        IHS_unstandardized, IHS_standardized=None):
        try:
            log.debug("inserting item...")
            item = IHSInfo.insert( population_set = population_set, 
                                locus               = locus, 
                                ones_freq           = float(ones_freq),
                                IHH_1               = float(IHH_1),
                                IHH_0               = float(IHH_0),
                                IHS_unstandardized  = float(IHS_unstandardized),
                                IHS_standardized    = float(IHS_standardized) if 
                                        IHS_standardized else None)
            return item
        except pw.IntegrityError:
            log.warning("iHS already in DB for locus {}, skipping.".format(
                                str(locus.variant_id)) )

    def add_xpehh_info(self, locus, population_set_a, population_set_b, 
        population_set_a_one_freq, IHH_a, population_set_b_one_freq, IHH_b, 
        XPEHH_unstandardized, XPEHH_standardized=None):

        item, created = XPEHHInfo.get_or_create( locus = locus,
                population_set_a          = population_set_a,
                population_set_b          = population_set_b,
                population_set_a_one_freq = float(population_set_a_one_freq),
                IHH_a                     = float(IHH_a),
                population_set_b_one_freq = float(population_set_b_one_freq),
                IHH_b                     = float(IHH_b),
                XPEHH_standardized        = float(XPEHH_standardized),
                XPEHH_unstandardized      = float(XPEHH_unstandardized) if 
                                                XPEHH_standardized else None)

        return item

    def get_popset_for_pop_names(self, population_name_list, 
                                                    superpopulation_name_list):
        '''
            Get the population set corresponding to the list of population names 
            specified in name_list.
        '''
        # the "<<" operator signifies an "IN" query
        pops = Population.select(Population.id).where(
                                Population.name << population_name_list)
        superpops = SuperPopulation.select(SuperPopulation.id).where(
                            SuperPopulation.name << superpopulation_name_list)

        pop_ids = []
        for pop in pops:
            pop_ids.append(int(pop.id))

        superpop_ids = []
        for superpop in superpops:
            superpop_ids.append(int(superpop.id))

        query = PopulationSet.select()

        popSets = {}
        for row in query:
            setPopulationIDs = []
            for pop in row.populations:
                setPopulationIDs.append(pop.id)

            setSuperpopulationIDs = []
            for superpop in row.super_populations:
                setSuperpopulationIDs.append(superpop.id)

            popSets[row.id] = { 
                                "pops":set(setPopulationIDs), 
                                "superpops":set(setSuperpopulationIDs)
                              }

        matchingPopSets = []
        for k,v in popSets.iteritems():
            if v["pops"]==set(pop_ids) and v["superpops"] == set(superpop_ids):
                matchingPopSets.append(k)


        if len(matchingPopSets):
            matchingPopSet = PopulationSet.get( 
                                        PopulationSet.id==matchingPopSets[0])

            return matchingPopSet
        else:
            return None

    def get_locus_info(self, chrom, pos_bp, locus_id=None):
        if not locus_id:
            return LocusInfo.get((LocusInfo.chrom == chrom) 
                                            & (LocusInfo.pos_bp == pos_bp))
        else:
            return LocusInfo.get((LocusInfo.variant_id == locus_id) 
                                            & (LocusInfo.chrom == chrom))

    def get_allele_info(self, locus, population_set):
        return AlleleInfo.get( (AlleleInfo.locus.id == locus.id) 
                        & (AlleleInfo.population_set.id == population_set.id) )

    def get_ehh_info(self, core_locus, population_set):
        return EHHInfo.select().where( (EHHInfo.core_locus.id == core_locus.id) 
                            & (EHHInfo.population_set.id == population_set.id) )

    def get_ihs_info(self, locus, population_set):
        return IHSInfo.get( (IHSInfo.locus.id == locus.id) 
                            & (IHSInfo.population_set.id == population_set.id) )        

    def get_xpehh_info(self, locus, population_set_a, population_set_b):
        return XPEHHInfo.get( (XPEHHInfo.locus.id == locus.id) & 
            (XPEHHInfo.population_set_a.id == population_set_a.id) & 
            (XPEHHInfo.population_set_b.id == population_set_b.id))

    def update_locus_info(self, chrom, variant_id, pos_bp, map_pos_cm, 
                            ref_allele, alt_allele, ancestral_call):
        if not locus_id:
            query = LocusInfo.update( chrom          = int(chrom),
                                      variant_id     = variant_id, 
                                      pos_bp         = int(pos_bp),
                                      map_pos_cm     = float(map_pos_cm),
                                      ref_allele     = ref_allele, 
                                      alt_allele     = alt_allele, 
                                      ancestral_call = ancestral_call 
                                      ).where((LocusInfo.chrom == chrom) 
                                        & (LocusInfo.pos_bp == pos_bp))
        else:
            locus = LocusInfo.update( variant_id     = variant_id, 
                                      pos_bp         = int(pos_bp),
                                      map_pos_cm     = float(map_pos_cm),
                                      ref_allele     = ref_allele, 
                                      alt_allele     = alt_allele, 
                                      ancestral_call = ancestral_call 
                                      ).where((LocusInfo.variant_id == locus_id) 
                                        & (LocusInfo.chrom == chrom))

        query.execute()

    def update_ehh_info(self, population_set, core_locus, pos_bp_delta, 
                            map_pos_cm_delta, EHH_ones, EHH_zeros, EHH):
        query = EHHInfo.update( population_set  = population_set, 
                                core_locus       = core_locus, 
                                pos_bp_delta     = int(pos_bp_delta),
                                map_pos_cm_delta = float(map_pos_cm_delta),
                                EHH_ones         = float(EHH_ones),
                                EHH_zeros        = float(EHH_zeros),
                                EHH              = float(EHH)
                               ).where( (EHHInfo.core_locus.id == core_locus.id) 
                             & (EHHInfo.population_set.id == population_set.id))

        query.execute()

    def update_ihs_info(self, population_set, locus, ones_freq, IHH_1, IHH_0, 
                        IHS_unstandardized, IHS_standardized=None):
        query = IHSInfo.update( population_set = population_set, 
                                locus               = locus, 
                                ones_freq           = float(ones_freq),
                                IHH_1               = float(IHH_1),
                                IHH_0               = float(IHH_0),
                                IHS_unstandardized  = float(IHS_unstandardized),
                                IHS_standardized    = float(IHS_standardized)
                                ).where( (IHSInfo.locus.id == locus.id) 
                             & (IHSInfo.population_set.id == population_set.id))

        query.execute()

    def update_xpehh_info(self, locus, population_set_a, population_set_b, 
                            population_set_a_one_freq, IHH_a, 
                            population_set_b_one_freq, IHH_b, 
                            XPEHH_unstandardized, XPEHH_standardized=None):
        query = XPEHHInfo.update( locus = locus,
                population_set_a          = population_set_a,
                population_set_b          = population_set_b,
                population_set_a_one_freq = float(population_set_a_one_freq),
                IHH_a                     = float(IHH_a),
                population_set_b_one_freq = float(population_set_b_one_freq),
                IHH_b                     = float(IHH_b),
                XPEHH_standardized        = float(XPEHH_standardized),
                XPEHH_unstandardized      = float(XPEHH_unstandardized)
                ).where( (XPEHHInfo.locus.id == locus.id) & 
            (XPEHHInfo.population_set_a.id == population_set_a.id) & 
            (XPEHHInfo.population_set_b.id == population_set_b.id))

        query.execute()

class CalculationReader(object):
    def __init__(self, jsonMetadataFile):
        self.jsonMetadataFile = os.path.abspath( 
                                        os.path.expanduser(jsonMetadataFile) )
        self.metadata = util.json_helpers.JSONHelper.read_data( 
                                                        self.jsonMetadataFile )

    def read_locus_info(self):
        if "tped_allele_metadata_file" in self.metadata:
            inFilePath = self.metadata["tped_allele_metadata_file"]
            if os.path.isfile( inFilePath ):
                with util.file.open_or_gzopen(inFilePath,"r") as inFile:
                    for line in csv.DictReader(inFile, delimiter="\t"):
                        yield line

    def read_ehh(self):
        if "ehh" in self.metadata:
            for locusDict in self.metadata["ehh"]:
                if "locus" in locusDict:
                    inFilePath = locusDict["ehh"]
                    if os.path.isfile( inFilePath ):
                        with open(inFilePath, "r") as inFile:
                            for line in csv.DictReader(inFile, 
                                            fieldnames=("physical_pos_bp_delta", 
                                            "genetic_pos_cm_delta", 
                                            "derived_1_ehh", 
                                            "ancestral_0_ehh", 
                                            "ehh"), 
                                            delimiter="\t"):
                                line["locus"] = locusDict["locus"]
                                yield line
                        

    def read_ihs(self):
        if "ihs" in self.metadata:
            inFilePath = self.metadata["ihs"]
            if os.path.isfile( inFilePath ):
                with open(inFilePath, "r") as inFile:
                    for line in csv.DictReader(inFile, 
                                        fieldnames=("locus_id",
                                        "physical_pos_bp",
                                        "ones_freq",
                                        "ihh1",
                                        "ihh0",
                                        "ihs_unstandardized"), 
                                        delimiter="\t"):
                        yield line

    def read_xpehh(self):
        if "xpehh" in self.metadata:
            for xpehh_scan in self.metadata["xpehh"]:
                inFilePath = xpehh_scan["xpehh"]
                if os.path.isfile( inFilePath ):
                    with open(inFilePath, "r") as inFile:
                        for line in csv.DictReader(inFile, delimiter="\t"):
                            yield line

class ScanStatStorer(CalculationReader):
    def __init__(self, metadata_json_file_path, db_path):
        # py2/3 compatible call of superclass constructor
        super(ScanStatStorer, self).__init__(metadata_json_file_path)

        self.metadata_json_file_path = metadata_json_file_path
        self.db_path = db_path

        self.dbm = DatabaseManager(db_path)
        self.store_population_set()

    def store_population_set(self):
        self.population_set = self.dbm.get_popset_for_pop_names(
            population_name_list=self.metadata["populations"], 
            superpopulation_name_list=self.metadata["super_populations"])
        if not self.population_set:
            self.population_set = self.dbm.add_population_set(
                population_list=self.metadata["populations"], 
                superpopulation_list=self.metadata["super_populations"])

    def store_locus_info(self):
        with self.dbm.db.atomic():
            for line in self.read_locus_info():
                self.dbm.add_locus_info( line["CHROM"],
                    line["VARIANT_ID"], 
                    line["POS_BP"],
                    line["MAP_POS_CM"],
                    line["REF_ALLELE"], 
                    line["ALT_ALLELE"], 
                    line["ANCESTRAL_CALL"] )

    def store_ihs(self):
        chromNum = self.metadata["chromosome_num"][0]

        for chunk in boltons.iterutils.chunked_iter(self.read_ihs(), 500):
            itemsToAdd = []
            locusIds = []
            locusIdsToIndices = {}
            for line in chunk:
                locusIds.append(line["locus_id"])

                itemsToAdd.append({"population_set" : self.population_set, 
                                "locus"               : line["locus_id"],
                                "ones_freq"           : float(line["ones_freq"]),
                                "IHH_1"               : float(line["ihh1"]),
                                "IHH_0"               : float(line["ihh0"]),
                                "IHS_unstandardized"  : float(line["ihs_unstandardized"]),
                                "IHS_standardized"    : None})

            # get list of loci matching IDs specified in this chunk
            with self.dbm.db.atomic():
                query = LocusInfo.select().where(
                    (LocusInfo.variant_id << locusIds) 
                    & (LocusInfo.chrom == chromNum))

            # build a dict mapping variant_id -> index
            for locusResult in query:
                #print locusResult.id, locusResult.variant_id
                locusIdsToIndices[locusResult.variant_id] = locusResult.id

            # set the ID of the locus in the list of items to insert
            # (faster than doing many selects to find and set the foreign keys)
            for item in itemsToAdd:
                item["locus"] = locusIdsToIndices[item["locus"]]

            with self.dbm.db.atomic():
                IHSInfo.insert_many(itemsToAdd).execute()

    def store_ehh(self):
        chromNum = self.metadata["chromosome_num"][0]
        
        with self.dbm.db.atomic():
            for line in self.read_ehh():
                locus = self.dbm.get_locus_info(chromNum, 0,
                                                        locus_id=line["locus"])

                self.dbm.add_ehh_info( self.population_set, 
                    locus, 
                    line["physical_pos_bp_delta"],
                    line["genetic_pos_cm_delta"],
                    line["derived_1_ehh"], 
                    line["ancestral_0_ehh"], 
                    line["ehh"] )

    def _store_pop_and_locus_for_other_xpehh_populations(self):
        if "xpehh" in self.metadata:
            for xpehh_scan in self.metadata["xpehh"]:
                metadata_file_for_pop2_tped = \
                 xpehh_scan["xpehh_pop_B_tped"].replace(
                                                    ".tped.gz",".metadata.json")
                cr = ScanStatStorer(metadata_file_for_pop2_tped, self.db_path)
                # population_set info is set by constructor of this class

                cr.store_locus_info()

                yield cr

    def store_xpehh(self):
        chromNum = self.metadata["chromosome_num"][0]

        for cr in self._store_pop_and_locus_for_other_xpehh_populations():
            for chunk in boltons.iterutils.chunked_iter(cr.read_xpehh(), 500):
                itemsToAdd = []
                locusIds   = []
                locusIdsToIndices = {}
                for line in chunk:
                    if line["pos"] not in locusIds:
                        locusIds.append(line["pos"])

                        itemsToAdd.append({"locus" : line["pos"],
                        "population_set_a"          : self.population_set,
                        "population_set_b"          : cr.population_set,
                        "population_set_a_one_freq" : line["p1"],
                        "IHH_a"                     : line["ihh1"],
                        "population_set_b_one_freq" : line["p2"],
                        "IHH_b"                     : line["ihh2"],
                        "XPEHH_unstandardized"      : line["xpehh"]})

                # get list of loci matching IDs specified in this chunk
                with self.dbm.db.atomic():
                    query = LocusInfo.select().where(
                        (LocusInfo.variant_id << locusIds) 
                        & (LocusInfo.chrom == chromNum))

                # build a dict mapping variant_id -> index
                for locusResult in query:
                    locusIdsToIndices[locusResult.variant_id] = locusResult.id

                # set the ID of the locus in the list of items to insert
                # (faster than doing selects to find and set the foreign keys)
                for item in itemsToAdd:
                    item["locus"] = locusIdsToIndices[item["locus"]]

                with self.dbm.db.atomic():
                    XPEHHInfo.insert_many(itemsToAdd).execute()

if __name__ == "__main__":
    pass



















