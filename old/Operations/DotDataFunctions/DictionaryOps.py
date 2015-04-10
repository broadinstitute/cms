
import operator

def MergeDicts( *dicts ):
    """Construct a merged dictionary from the given dicts.
    If two dicts define the same key, the key from the dict later in the list is chosen."""
    return dict( reduce( operator.add, map( dict.items, dicts ) ) )

def RestrictDict( aDict, restrictSet ):
    """Return a dict which has the mappings from the original dict only for keys in the given set"""
    return dict( [ ( k, aDict[k] ) for k in frozenset( restrictSet ) & frozenset( aDict.keys() ) ] )
