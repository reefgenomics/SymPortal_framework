class Error(Exception):
   """Base class for other exceptions"""
   pass

class InsufficientSequencesInAlignment(Error):
    """Raised when there are insufficient sequences for making a tree"""
    pass