class CMSException(Exception):
    '''base class for CMS exceptions. except for built-ins, all should be descendants of this exception so that an except block for this exception will catch all others'''
    pass

class CMSInputError(CMSException):
    pass
