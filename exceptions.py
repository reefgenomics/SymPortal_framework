class Error(Exception):
   """Base class for other exceptions"""
   pass

class InsufficientSequencesInAlignment(Error):
    """Raised when there are insufficient sequences for making a tree"""
    pass

class NoDataSetSampleSequencePMObjects(Error):
    """Raised when there are no DataSetSampleSequencePM objects associated
    to a sample. This will only be raised in the output context if we are mistakenly
    trying to generate a pre med sequence output table but there DataSetSampleSequencePM
    objects were not generated during data loading"""
    pass