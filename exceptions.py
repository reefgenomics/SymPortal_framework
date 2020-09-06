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

class DistanceTypeNotIdentifiedError(Error):
    """Raised when looking for the distance type i.e. unifrac or braycurtis in the path of the PCoA csv.
    Will also be raised if the transofrmation type cannot be inferred i.e. sqrt or no_sqrt"""
    pass

class EigenValsTooSmallError(Error):
    """Raised when the eigen values of a PCoA sum to 0. This happens because
    the scipy code converts very small eigen values to 0. As such if all eigen values are small
    we end up with TrueDivide errors when trying to calculate variance."""
    pass