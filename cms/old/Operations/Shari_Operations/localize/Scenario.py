
from Operations.MiscUtil import namedtuple, ApplyToResult
from Operations.Shari_Operations.localize.PopConsts import AllAges, AllPops, AllFreqs
import re

class Scenario(object):
    """Represents one simulation scenario."""

    @staticmethod
    def fromString( s, mutAge = 10 ):
        """Return a scenario based on its string representation"""
        if isinstance( s, Scenario ): return s
        if s == 'neutral': return NeutralScenario()

        m = re.match( '(\d+)ky/sel(\d+)_(\d)', s )
        if m: return SelectionScenario( mutAge = int( m.group( 1 ) ), mutFreq = int( m.group( 2 ) ), mutPop = int( m.group( 3 ) ) )
        
        m = re.match( 'sel(\d+)_(\d)', s )
        if m: return SelectionScenario( mutAge = mutAge, mutFreq = int( m.group( 1 ) ), mutPop = int( m.group( 2 ) ) )
        raise ValueError( 'invalid scenario - ' + s )

    fromStr = fromString
            
    def isNeutral( self ): return self.is_neutral()

class NeutralScenario(Scenario, namedtuple('NeutralScenario', '')):
    """Represents a simulation of neutral evolution ( no selection )"""
    def is_neutral( self ): return True
    def scenName( self ):
        """Return a string representing the scenario"""
        return 'neutral'
    def scenDir( self ):
        """Return a string representing a subdir in which data for this scenario is stored."""
        return 'neutral'

    def __bool__(self): return True

    def __str__( self ):
        """Return an informal string representation"""
        return self.scenDir()

    def __repr__( self ):
        """Return a formal string representation"""
        return 'NeutralScenario()'
        
    

class SelectionScenario(Scenario, namedtuple('SelectionScenario', 'mutAge mutPop mutFreq')):
    """Represents a simulation of evolution where there was a selective event.

    Fields:

      mutAge - how long ago the good mutation arose ( in kiloyears )
      mutPop - in which population the good mutation arose.  this is an integer.
        currently, the integer is 1,4,5 with the correspondence {1:'European',4:'EastAsian',5:'WestAfrican'}.
      mutFreq - frequency of this mutation _today_ (what fraction of the population
         in which the mutation arose, has this mutation today).
         Note that this is only the average frequency of replicas simulated with this scenario!
         In each simulated replica, the exact present-day frequency of the selected allele in the selected population
         will differ, though hopefully the values across all replicas will cluster around the scenario's
         specified frequency.
      
    """
    def is_neutral( self ): return False
    def scenName( self ):
        """Return a string representing the scenario (does not include age)."""
        return 'sel%d_%d' % ( self.mutFreq, self.mutPop )
    def scenDir( self ):
        """Return a string representing a subdir in which data for this scenario is stored."""
        return '%sky/%s' % ( self.mutAge, self.scenName() )
    def __str__( self ):
        """Return an informal string representation"""
        return self.scenDir()

    def __repr__( self ):
        """Return a formal string representation"""
        return 'SelectionScenario( mutAge = %d, mutPop = %d, mutFreq = %d )' % ( self.mutAge, self.mutPop, self.mutFreq )

@ApplyToResult(tuple)    
def GetScenarios( mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs, includeNeutral = True ):
    """Yield the neutral scenario plus the full set of
    selection scenarios for all combinations of (mutAge, mutPop, mutFreq).  The neutral scenario is always
    yielded first (some code relies on this)."""
    if includeNeutral: yield NeutralScenario()
    for mutAge in mutAges:
        for mutPop in mutPops:
            for mutFreq in mutFreqs:
                yield SelectionScenario( mutAge, mutPop, mutFreq )

@ApplyToResult(tuple)
def GetSelectionScenarios( mutAges = AllAges, mutPops = AllPops, mutFreqs = AllFreqs ):
    """Yield the selection scenarios for all combinations of (mutAge, mutPop, mutFreq)"""
    return GetScenarios( mutAges, mutPops, mutFreqs, includeNeutral = False )

ScenarioSetP = namedtuple( 'ScenarioSet', 'mutAges mutPops mutFreqs scenarios' )
class ScenarioSet( ScenarioSetP ):
    """Represents a set of scenarios that includes the neutral scenario,
    as well a selection scenarios for all combinations of the given
    (mutAges mutPops mutFreqs).
    """

    def __new__( cls, mutAges, mutPops, mutFreqs ):
        return ScenarioSetP.__new__( cls, mutAges, mutPops, mutFreqs,
                                     scenarios = tuple( GetScenarios( mutAges, mutPops, mutFreqs ) ) )
    
    

