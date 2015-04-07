from numpy import longdouble, ulonglong, sqrt

class RunningStat(object):
    """Class for computing running mean and stddev.  The class uses only a constant amount of space.

    Usage:

    Call add() to add values.  Call var(), std() and mean() to get the
    current running statistics.  Note that the values are returned as
    type numpy.longdouble .

    """

    def __init__( self ):

        self.s = longdouble( 0 )
        self.m = longdouble( 0 )
        self.last_m = longdouble( 0 )
        self.n = ulonglong( 0 )
        self.is_started = False

    def add( self, x ):
        """Add a value, and update the running stats"""
        self.n += 1
        x = longdouble( x )
        if not self.is_started:
            self.m = x
            self.s = longdouble( 0 )
            self.is_started = True
        else:
            self.last_m = self.m
            self.m += ( x - self.m ) / longdouble( self.n )
            self.s += ( x - self.last_m ) * ( x - self.m )

    def var( self ):
        """Return the running variance"""
        return self.s / longdouble( self.n )

    def std( self ):
        """Return the running standard deviation"""
        return sqrt( self.var() )

    def mean( self ):
        """Return the running mean"""
        return self.m
            


