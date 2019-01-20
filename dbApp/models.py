from django.db import models
from django.db.models import Max
import json

# python3 manage.py graph_models -a -g -o my_project.svg
# You can visualise these models using the following commands in the terminal
# http://cheng.logdown.com/posts/2015/07/07/visualize-the-relationships-of-django-models-using-django-extensions
# sudo apt install graphviz I had to install that in the terminal as well (normal terminal) i.e. not with python3 etc.


class SymportalFramework(models.Model):
    latest_reference_fasta = models.CharField(max_length=30, default='symClade_2_2.fa')
    # this will be fed to the screen_sub_e_value_sequences method when making the next database iteration
    next_reference_fasta_iteration = models.IntegerField(default=1)
    # a sub e value sequence is one that matched one of the sequence in the reference fasta
    # but was thrown out because the evalue support for the match was not high enough
    # during each data_set, these sequences are kept track of.
    # If the sequences are found in a number of samples, then they are considered not to be sequencing artefacts
    # this number is the required_sub_e_value_seq_support below.
    required_sub_e_value_seq_support_samples = models.IntegerField(default=3)
    # sequences that are found in the above number of samples are then screened against the blast nt database
    # for their top 10 matches. If a number of matches are Symbiodinium then they are considered Symbiodinium
    # derived sequences and they are added to the reference fasta
    required_sub_e_value_seq_support_blast_symbiodinium = models.IntegerField(null=2)


class DataSet(models.Model):
    name = models.CharField(max_length=60, default='something')
    reference_fasta_database_used = models.CharField(max_length=60, default='None')
    submittingUser = models.CharField(max_length=100, default='no_user_defined')
    submitting_user_email = models.CharField(max_length=100, default='no_email_defined')
    # submittedDataRaw = models.FileField()
    workingDirectory = models.CharField(max_length=300, default='None')

    # This is true when all processing is complete, including MED analyses
    dataProcessed = models.BooleanField(default=False)
    # This is true when all sample processing is complete, i.e. MED may not necessarily be complete yet.
    initialDataProcessed = models.BooleanField(default=False)
    # This will be true if there is a python process currently processing this data_set
    # NB as the program may crash at any time this may get stuck to True
    # To prevent this we will have to set all of these to false on start up of SymPortal
    # On startup all dataSubs have their currentlyBeingProcessed to False using the apps.py and __init__.py scipts:
    # http://stackoverflow.com/questions/28896627/executing-code-on-startup-in-django-1-7
    currentlyBeingProcessed = models.BooleanField(default=False)
    timeStamp = models.CharField(max_length=100, default='None')
    # http://stackoverflow.com/questions/27220480/django-datetime-migration-error
    # We can programatically set whenSubmitted by '= datetime.now()'
    # The Datetimes were causing a lot of problems when doing the dumpdata and loaddata. I think it was becuase of the
    # for mat of the date times, specifically because they had a time zone applied.
    # I think we can likely fix this by specifying formats for the datetimes in the settings file
    # We can then make sure thatn when we enter a datetime using the python that we use the correct format.
    # https://docs.djangoproject.com/en/1.9/ref/settings/#datetime-input-formats

    def __str__(self):
        return self.name

    # http://stackoverflow.com/questions/11210396/django-set-field-values-dynamically-based-on-other-fields
    # http://stackoverflow.com/questions/4380879/django-model-field-default-based-off-another-field-in-same-model
    def save(self, *args, **kwargs):
        # This used to have content for percent complete to be calculated but I have gotten rid of this
        # I have left this structure here in case we want to use it again.
        # Content would go here
        super().save(*args, **kwargs)


class DataSetSample(models.Model):

    dataSubmissionFrom = models.ForeignKey(DataSet, on_delete=models.CASCADE, null=True)
    name = models.CharField(max_length=200, default='None')
    # This is the absolute number of sequences after make.contigs
    initialTotSeqNum = models.IntegerField(default=0)

    # store the aboslute number of sequences after sequencing QC at this stage
    post_seq_qc_absolute_num_seqs = models.IntegerField(default=0)
    # This is the unique number of sequences after the sequencing QC
    initialUniqueSeqNum = models.IntegerField(default=0)

    # Absolute number of sequences after sequencing QC and screening for Symbiodinium (i.e. Symbiodinium only)
    finalTotSeqNum = models.IntegerField(default=0)
    # Same as above but the number of unique seqs
    finalUniqueSeqNum = models.IntegerField(default=0)

    # store the abosolute number of sequenes that were not considered Symbiodinium
    non_sym_absolute_num_seqs = models.IntegerField(default=0)
    # This is the number of unique sequences that were not considered Symbiodinium
    nonSymSeqsNum = models.IntegerField(default=0)

    # store the abosulte number of sequences that were lost during the size selection
    size_violation_absolute = models.IntegerField(default=0)
    # store the unique number of sequences that were lost during the size screening
    size_violation_unique = models.IntegerField(default=0)

    # store the number of absolute sequences remaining after MED
    post_med_absolute = models.IntegerField(default=0)
    # store the number of unique sequences remaining after MED (nodes)
    post_med_unique = models.IntegerField(default=0)

    # Whether processed upto cladal separation but not MED yet
    initialProcessingComplete = models.BooleanField(default=False)
    # Whether processed as above but also including MED
    finalProcessingComplete = models.BooleanField(default=False)
    errorInProcessing = models.BooleanField(default=False)
    errorReason = models.CharField(max_length=100, default='noError')
    cladalSeqTotals = models.CharField(max_length=5000, null=True)

    # Meta data for the sample
    sample_type = models.CharField(max_length=50, default='NoData')
    host_phylum = models.CharField(max_length=50, default='NoData')
    host_class = models.CharField(max_length=50, default='NoData')
    host_order = models.CharField(max_length=50, default='NoData')
    host_family = models.CharField(max_length=50, default='NoData')
    host_genus = models.CharField(max_length=50, default='NoData')
    host_species = models.CharField(max_length=50, default='NoData')
    collection_latitude = models.DecimalField(max_digits=11, decimal_places=8, default=999.99999999)
    collection_longitude = models.DecimalField(max_digits=11, decimal_places=8, default=999.99999999)
    # do not use the django date field as this causes problems when trying to dump and load the database
    collection_date = models.CharField(max_length=40, default='NoData')
    # store a string rather than a number as this may be given as a range e.g. 6 - 12
    collection_depth = models.CharField(max_length=40, default='NoData')

    def __str__(self):
        return self.name


class DataAnalysis(models.Model):
    # This will be a jsoned list of uids of the dataSubmissions that are included in this analysis
    listOfDataSubmissions = models.CharField(max_length=500, null=True)
    withinCladeCutOff = models.FloatField(default=0.04)
    typeSupport = models.FloatField(default=0.01)
    cladeCollectionPopulationComplete = models.BooleanField(default=False)
    analysisTypesDefined = models.BooleanField(default=False)
    initialTypeDiscoComplete = models.BooleanField(default=False)
    analysisTypesAssigned = models.BooleanField(default=False)
    analysisTypesCollapsed = models.BooleanField(default=False)
    refSeqsNamed = models.BooleanField(default=False)
    speciesAssociated = models.BooleanField(default=False)
    name = models.CharField(max_length=100, null=True)
    description = models.CharField(max_length=5000, null=True)
    dbVersionAnalysis = models.BooleanField(default=False)
    dbVersionToRunAgainstID = models.IntegerField(null=True)
    dSIDToAnalyse = models.IntegerField(null=True)
    timeStamp = models.CharField(max_length=100, default='None')
    submittingUser = models.CharField(max_length=100, default='no_user_defined')
    submitting_user_email = models.CharField(max_length=100, default='no_email_defined')

    def get_clade_collections(self):
        list_of_uids = [int(x) for x in self.listOfDataSubmissions.split(',')]
        clade_collections = CladeCollection.objects.filter(dataSetSampleFrom__dataSubmissionFrom__in=list_of_uids)
        return clade_collections


class CladeCollection(models.Model):
    dataSetSampleFrom = models.ForeignKey(DataSetSample, on_delete=models.CASCADE, null=True)
    clade = models.CharField(max_length=1)
    # the method below to get the footprint of the clade_collection_object is incredibly slow.
    # Hopefully a much faster way will be to store the ID's of the refseqs that make up the footprint
    # I will therefore store the reference uids in the field below
    footPrint = models.CharField(max_length=100000, default=True)

    # maj() should return the analysedSampleSequence associated with the
    # cladeCollection in Q with the largest abundance
    def maj(self):
        max_abundance = DataSetSampleSequence.objects.filter(
            cladeCollectionTwoFoundIn=self).aggregate(Max('abundance'))['abundance__max']
        # A .get() could fail if there are two sequences with an equal abundance so I
        # will try the get first and then failing that I will do a .filter() and return the first element
        try:
            return DataSetSampleSequence.objects.get(cladeCollectionTwoFoundIn=self, abundance=max_abundance)

        except:
            return DataSetSampleSequence.objects.filter(cladeCollectionTwoFoundIn=self, abundance=max_abundance)[0]

    # range will be how many of the seqs we want to return
    def ordered_list_of_reference_sequences(self, num_to_return):
        if num_to_return:
            dsss = list(DataSetSampleSequence.objects.filter(cladeCollectionTwoFoundIn=self).order_by('-abundance'))
            ordered_reference_sequences = [ds.referenceSequenceOf for ds in dsss]
        else:
            # Try is for incase there are not 10 seqs in the clade_collection_object
            try:
                dsss = list(DataSetSampleSequence.objects.filter(
                    cladeCollectionTwoFoundIn=self).order_by('-abundance')[:10])
                ordered_reference_sequences = [ds.referenceSequenceOf for ds in dsss]
            except:
                dsss = list(DataSetSampleSequence.objects.filter(
                    cladeCollectionTwoFoundIn=self).order_by('-abundance'))
                ordered_reference_sequences = [ds.referenceSequenceOf for ds in dsss]
        return ordered_reference_sequences

    # This will return the foot print of the analysedSampleSequences that are found above the given percentage cutoff
    def cutoff_footprint(self, cutoff):
        # get total seqs in cladeCollection
        total = 0

        for dsss in DataSetSampleSequence.objects.filter(cladeCollectionTwoFoundIn=self):
            total += dsss.abundance
        sequence_number_cutoff = cutoff * total
        frset = set([])
        for dsss in DataSetSampleSequence.objects.filter(cladeCollectionTwoFoundIn=self):
            if dsss.abundance > sequence_number_cutoff:
                frset.add(dsss.referenceSequenceOf)
        return frozenset(frset)

    # footprint should return a frozen of referenceSequenceObjects associated
    # with the analysedSampleSequences associated with the cladeCollectionInQ
    def footprint(self):
        return frozenset(
            refseq for refseq in ReferenceSequence.objects.filter(id__in=[int(a) for a in self.footPrint.split(',')]))

    def __str__(self):
        return self.dataSetSampleFrom.name


class AnalysisGroup(models.Model):
    dataAnalysisFrom = models.ForeignKey(DataAnalysis, on_delete=models.CASCADE, null=True)
    name = models.CharField(max_length=100, null=True)

    # We will generate a name for the analysis_group at some point later according to the analysisTypes it contains


class AnalysisType(models.Model):
    dataAnalysisFrom = models.ForeignKey(DataAnalysis, on_delete=models.CASCADE, null=True)
    analysisGroupOf = models.ForeignKey(AnalysisGroup, on_delete=models.SET_NULL, null=True)
    # This should be a frozen set of referenceSequences
    # As this is not possible in a django field or jsonable
    # I will instead make it a string of reference sequence uids
    # I will make set and get methods for the footprint
    # This will be a list of refSeqs that make up the footprint in order of their abundance when type first defined
    # This is a commar separated string of the uids of the ref seqs that define the type
    orderedFootprintList = models.CharField(max_length=200, null=True)
    # Same for listOfMajs
    # set() of refseqs that are Majs in each of the CCs this type was initially identified in.
    # Note that this is therefore in no particular order
    MajRefSeqSet = models.CharField(max_length=40, null=True)
    # Same for this
    # The list of cladeCollections in which the type was defined from
    listOfCladeCollectionsFoundInInitially = models.CharField(max_length=5000, null=True)
    # The list of cladeCollections that the type was associated with after one iteration of assigningTypes
    listOfCladeCollections = models.CharField(max_length=100000, null=True)
    # This is a 2D list, a list for each clade collection in order of the listofCladeCollections
    # Within each list the absolute abundances of the defining seqs in order of orderedFootprintList
    footprintSeqAbundances = models.CharField(max_length=100000, null=True)
    # Same as above but the relative abundance of the
    # seq in Q as a function of all of the sequences in the cladeCollection
    footprintSeqRatios = models.CharField(max_length=100000, null=True)
    clade = models.CharField(max_length=1)
    coDom = models.BooleanField(default=False)

    name = models.CharField(max_length=1000, null=True)
    # [(max,min) for refseq in ordered footprint list]
    maxMinRatios = models.CharField(max_length=100000, null=True)
    # The list of speceis that this type is associated with
    species = models.CharField(max_length=200, null=True)

    # this list will keep track of which of the defining intras of this type are 'unlocked' i.e. at least one
    # of the instances of that intra were found at <5%. We will use this list to the 'artefact type creation'
    # This artefactIntras will hold a char string of comma separated ints that represent the id's of the
    # refseqs that are unlocked
    artefactIntras = models.CharField(max_length=5000, default='')

    # This will only be set to true when running an analysis against a database version.
    isLockedType = models.BooleanField(default=False)

    # This method will populate the following attributes
    #  self.listOfCladeCollectionsFoundInInitially
    # self.orderedFootprintList
    # self.name
    # self.footprintSeqAbundances
    # Very important that the listOfCC be an acutal list rather than a query set as the query set may not maintain order
    # We need order to be maintained when we come to do the multimodal detection
    # I have checked and the listOfCC is indeed a list rather than a queryset
    # This creates the initial seqabundance information, maxmin ratios and seqratios for use
    # in the multimodal detection and the ratio checking for type association.
    # Once we have been through the initial round of type association we will
    # be using the updateTypeAttribute method in place of this so that we will start working from the
    # listofCladeCaoolections rather than the listOfCladeCollectionsFoundInInitially
    # TODO because I am going to move away from this silly ratio business and work on relative abundance
    # of defining intras.
    # These will simply be the relative proportion of any given intra (including the maj) with respect to the
    # total number of sequences the intras of the type make up.

    def init_type_attributes(self, list_of_clade_collections, footprintlistofrefseqs):
        # Becuase this can be run twice on a type due to the checking for within clade cutoff methods
        # We need to reset the artefactIntras
        self.artefactIntras = ''

        self.listOfCladeCollectionsFoundInInitially = ','.join([str(cc.id) for cc in list(list_of_clade_collections)])

        footprintlist = list(footprintlistofrefseqs)
        count_list = [[] for _ in list_of_clade_collections]

        # Here get the abundance of each of the defining intras in each of the samples of this type
        # Then we can sort it and make a name

        # here we are going to introduce a new concept.
        # it is to do within minimising the effect of the 3% withinCladeCuttoff.
        # An example case is maybe easiest.
        '''
        +----------------------+-------+-------+-------+-------+-------+-------+
        |         Type         |  D1   |  D4   |  D6   |  D17  | D2.2  |  D2   |
        +======================+=======+=======+=======+=======+=======+=======+
        | D4/D1/D6-D2-D2.2     | 0.335 | 0.320 | 0.182 | 0.017 | 0.072 | 0.071 |
        +----------------------+-------+-------+-------+-------+-------+-------+
        | D1/D4/D17-D6-D2.2-D2 | 0.394 | 0.221 | 0.107 | 0.178 | 0.061 | 0.036 |
        +----------------------+-------+-------+-------+-------+-------+-------+
        '''
        # The two above types were created in type discovery. The difference being a lack of D17. The types are roughly
        # carried through to the end. They are not collapsed due to the much higher average of D17 in the second group
        # However, it is not to say that the CCs assigned the first type don't contain D17 infact they do, only at
        # below the 3% cutoff that has been introduced by SP. However, there are CCs in the first second type that have
        # D17 at very low abundances, like 3 or 4 %. As such this difference in these types is representative of an
        # artefact rather than biology. To mitigate this effect when a profile is created, for each of its intras, in
        # each of its CCs, we will see if the intra is found at <= 0.05 retlative abundance i.e. 5%.
        # If at least one of the intras in one of the clade_collection_object in the type is found, then we
        # will consider this profile
        # to be at danger of being effected by the SP cut off artefact effect.
        # In this case, we will lower the lower boundary to 0.5%. From empirical data this seems to include most cases
        # there are some CCs that have hyper low amounds of it but this not turely representative of the sequence.
        # In this way we should mitigate the issues we are seeing.
        # I think the easiest way to implement this is here. Where we have an intra in a type as discussed above
        # we simply set the minimum value to 0.005 so that this value is automatically used.

        # I don't think that it is necessary to predetermine whether a type is CoDom or not
        # we do all of the work that is required for that here anyway.
        # So instead I am going to workout coDom in the type init here.
        # work out the maj seqs here
        set_of_majority_reference_sequences = set()
        for clade_collection_object in list_of_clade_collections:
            max_abund = 0
            maj_reference_sequence = None
            list_of_sequences_from_clade_collection = DataSetSampleSequence.objects.filter(
                cladeCollectionTwoFoundIn=clade_collection_object)
            for ref_seq in footprintlistofrefseqs:
                sample_sequence_of_reference_sequence = list_of_sequences_from_clade_collection.get(
                    referenceSequenceOf=ref_seq)
                abundance_of_reference_sequence_in_question = sample_sequence_of_reference_sequence.abundance
                if abundance_of_reference_sequence_in_question > max_abund:
                    maj_reference_sequence = sample_sequence_of_reference_sequence.referenceSequenceOf
                    max_abund = sample_sequence_of_reference_sequence.abundance
                count_list[
                    list_of_clade_collections.index(
                        clade_collection_object)].append(abundance_of_reference_sequence_in_question)
            set_of_majority_reference_sequences.add(maj_reference_sequence)

        # if set_of_majority_reference_sequences > 1 then this is coDom
        if len(set_of_majority_reference_sequences) > 1:
            self.coDom = True
        else:
            self.coDom = False
        self.set_maj_ref_seq_set(set_of_majority_reference_sequences)
        # At this point we have the abundance of each sequence for each instance of the type
        # we have 2d array where list are cladeCollections, and items in each list represent sequence
        # Now create counting list to count the abundance of each seq in the footprints refseqs
        abundance_list = [[refseq, 0] for refseq in footprintlist]
        for i in range(len(count_list)):  # for each cladeCollection with footprint
            for j in range(len(count_list[0])):  # for each refseq in footprint
                abundance_list[j][1] += count_list[i][j]

        # order the list
        ordered_footprint_list = [a[0] for a in sorted(abundance_list, key=lambda x: x[1], reverse=True)]
        # This now gives us the abundance order of the refseqs
        # TODO 08/12/17 this is currently working with absolute abundances instead of relative abundances
        # Which strikes me as incorrect. BUT I don't want to start messing with this now.
        # just so long as we work with the relative abundances as funtion of the
        # type seuences in the clade_collection_object

        self.orderedFootprintList = ','.join([str(refseq.id) for refseq in ordered_footprint_list])

        if not self.isLockedType:
            self.name = self.generate_name(ordered_footprint_list)

        # The below code re orders the countlist 2D list so that the columns are in order of refseq abundance
        count_list_ordered = []
        for row in count_list:
            count_list_ordered.append([])  # new empty list for each row
            for element in ordered_footprint_list:  # for each refseq in order of abundance
                count_list_ordered[-1].append(row[footprintlist.index(element)])

        # recover the memory from countList
        del count_list
        # Here we have the 2d list of raw abundances of seqs for every clade collection (row) for every seq (column)
        # Now we can store it in a char field by json.dumpsing it.
        self.footprintSeqAbundances = json.dumps(count_list_ordered)

        # These now represent relative abundances rather than ratios

        self.generate_max_min_ratios(count_list_ordered, ordered_footprint_list)
        self.generate_ratio_list()
        self.save()

    # This is similar to initTypeAttributes but it works on the listOfCladeCollections
    # It sets self.orderedFootprintList, self.name, self.footprintSeqAbundances, self.genreateMaxMinRatios,
    # self.coDom, self.MajRefSeqSet
    # the self.maxminratios will be kept uptodate through the assigning types to clade collections process but
    # This will be called everytime we finish an iteration of the assigning types. after the first iteration for
    # this will convert from working on the initial clade collection list to the listOfCladeCollections
    # It will be important to call this before running the miltimodal function as the counts and ratios will
    # need to be up to date.
    # The listOfCladeCollections should be kept up to date as we go through the association of types.

    def update_type_attributes(self):
        # Important to use this get() format rather than all().filter() format as
        # we need to maintain the order of the cladeCollection list
        # NB this is now updating from listoflcadecollections rather than from listoflcadecollectionsfoundininitially

        list_of_clade_collections = [
            CladeCollection.objects.get(id=ID) for ID in [int(x) for x in self.listOfCladeCollections.split(',')]]

        ordered_list_of_reference_sequences_in_type = [
            ReferenceSequence.objects.get(id=ID) for ID in [int(x) for x in self.orderedFootprintList.split(',')]]

        count_list = [[] for _ in list_of_clade_collections]

        for clade_collection_object in list_of_clade_collections:
            list_of_sequences_from_clade_collection = DataSetSampleSequence.objects.filter(
                cladeCollectionTwoFoundIn=clade_collection_object)
            for ref_seq in ordered_list_of_reference_sequences_in_type:
                sample_sequence_of_reference_sequence = list_of_sequences_from_clade_collection.get(
                    referenceSequenceOf=ref_seq)
                count_list[
                    list_of_clade_collections.index(
                        clade_collection_object)].append(sample_sequence_of_reference_sequence.abundance)

        # Here I will update/initiate self.MajRefSeqSet and self.coDom using the countList info
        majority_reference_sequence_set = set()
        for CClist in count_list:
            max_value = max(CClist)

            index = CClist.index(max_value)
            majority_reference_sequence_set.add(ordered_list_of_reference_sequences_in_type[index])
        self.MajRefSeqSet = ','.join([str(a.id) for a in list(majority_reference_sequence_set)])
        if len(majority_reference_sequence_set) > 1:
            self.coDom = True
        else:
            self.coDom = False
        self.save()
        # At this point we have the abundance of each sequence for each instance of the type
        # we have 2d array where list are samples and items in list represent sequence
        # Now create counting list to count the abundance of each seq in the footprints refseqs
        abundance_list = [[refseq, 0] for refseq in ordered_list_of_reference_sequences_in_type]
        for i in range(len(count_list)):  # for each cladeCollection with footprint
            for j in range(len(count_list[0])):  # for each refseq in footprint
                abundance_list[j][1] += count_list[i][j]

        ordered_footprint_list = [a[0] for a in sorted(abundance_list, key=lambda x: x[1], reverse=True)]
        # Because this is also used when creating a brand new type in the bimodal distribution splitting
        # We should always create a new name
        # if self.orderedFootprintList != ','.join([str(refseq.id) for refseq in orderedFootprintList]):
        # Then the order of the footprint has changed and the name needs to be redone
        # Else it doesn't need to be redone
        self.orderedFootprintList = ','.join([str(refseq.id) for refseq in ordered_footprint_list])

        if not self.isLockedType:
            self.name = self.generate_name(ordered_footprint_list)
        # If the order of the footPrint hasn't changed then the count list is
        # already in the right order and we don't need to do the code below.
        # If the order has changed then we do need to redo the code
        # The below code re orders the countlist 2D list os that the columns are in order of refseq abundance
        count_list_ordered = []
        for row in count_list:
            count_list_ordered.append([])  # new empty list for each row
            for element in ordered_footprint_list:  # for each refseq in order of abundance
                count_list_ordered[-1].append(row[ordered_list_of_reference_sequences_in_type.index(element)])
        # recover the memory from countList
        del count_list

        # Here we have the 2d list of raw abundances of seqs for every clade collection (row) for every seq (column)
        # Now we can store it in a char field by json.dumpsing it.

        self.footprintSeqAbundances = json.dumps(count_list_ordered)
        # These now represent relative abundances rather than ratios
        self.generate_max_min_ratios(count_list_ordered, ordered_footprint_list)
        self.generate_ratio_list()

    def add_clade_collection_to_clade_collection_list(self, clade_collection_object):
        if self.listOfCladeCollections:  # If empty then no need to add comma
            self.listOfCladeCollections += ',{}'.format(clade_collection_object.id)
        else:
            self.listOfCladeCollections = str(clade_collection_object.id)

    def add_clade_collection_list_to_clade_collection_list(self, clade_collection_uid_list):
        if self.listOfCladeCollections:  # If empty then no need to add commar
            self.listOfCladeCollections += ',{}'.format(','.join([str(a) for a in clade_collection_uid_list]))
        else:
            self.listOfCladeCollections = ','.join([str(a) for a in clade_collection_uid_list])

    def remove_clade_collection_list_from_initial_clade_collection_list(self, clade_collection_uid_list):
        list_of_uids_as_string = self.listOfCladeCollectionsFoundInInitially.split(',')
        new_list = [uid for uid in list_of_uids_as_string if uid not in clade_collection_uid_list]
        self.listOfCladeCollectionsFoundInInitially = ','.join(new_list)
        self.save()

    def remove_clade_collection_from_clade_collection_list(self, clade_collection_object):
        list_of_uids_as_string = self.listOfCladeCollections.split(',')
        list_of_uids_as_string.remove(str(clade_collection_object.id))
        self.listOfCladeCollections = ','.join(list_of_uids_as_string)
        self.save()

    # This is the relative abundance of each sequence in the type as a function of the total sequence found
    # in the clade Collection THAT ARE ALSO in the type. SO proportion of type sequences in cladecollection

    def generate_ratio_list(self):
        new_ratio_list = []
        footprint_seq_abundances_loaded = json.loads(self.footprintSeqAbundances)
        # For each clade_collection_object the type is associated with
        for i in range(len(footprint_seq_abundances_loaded)):
            new_ratio_list.append([])
            denominator = sum(footprint_seq_abundances_loaded[i])
            for j in range(len(footprint_seq_abundances_loaded[i])):  # For each intra in the type
                ratio = float("%.3f" % (footprint_seq_abundances_loaded[i][j]/denominator))
                new_ratio_list[-1].append(ratio)
        self.footprintSeqRatios = json.dumps(new_ratio_list)

    def get_ratio_list(self):
        return json.loads(self.footprintSeqRatios)

    def generate_name(self, ordered_list_of_reference_sequences_in_footprint=None, accession=None):
        # If self.isLockedType then we skip the naming.
        if not self.isLockedType:

            # this is interesting, sometimes one of the intras that is not a coDom is more abundant than one
            # of the coDom intras to solve this problem we should add the coDoms intas to the name first and then
            # go in the remaining order of the orderedListOfRefSeqsInFootprint

            if accession:
                # Here we output the accession name for the type in Q
                # This checks to see if there are accession numbers available for each of the refseqs in the type
                # If there are accessions available we use these, in the same / and - format as the usual name
                # Where there aren't accessions available for the ref seqs we use the ref seq's ID.
                orderedlistofrefseqsinfootprint = self.get_ordered_footprint_list()
                if self.coDom:
                    list_of_majority_reference_sequences = ReferenceSequence.objects.filter(
                        id__in=[int(x) for x in self.MajRefSeqSet.split(',')])
                    # Start the name with the coDom intras in order of abundance.
                    # Then append the noncoDom intras in order of abundance
                    ordered_list_of_majority_reference_sequences = [
                        ref_seq for ref_seq in orderedlistofrefseqsinfootprint if
                        ref_seq in list_of_majority_reference_sequences]

                    name = '/'.join(
                        [ref_seq.accession if ref_seq.accession is not None
                         else str(ref_seq.id) for ref_seq in ordered_list_of_majority_reference_sequences])
                    list_of_remaining_reference_sequences = [
                        ref_seq for ref_seq in orderedlistofrefseqsinfootprint if
                        ref_seq not in list_of_majority_reference_sequences]
                    if list_of_remaining_reference_sequences:
                        name += '-{}'.format('-'.join(
                            [ref_seq.accession if ref_seq.accession is not None else str(ref_seq.id) for ref_seq in
                             list_of_remaining_reference_sequences]))
                    return name
                else:
                    new_list = []
                    for ref_seq in orderedlistofrefseqsinfootprint:
                        if ref_seq.accession and ref_seq.accession is not None:
                            new_list.append(ref_seq.accession)
                        else:
                            new_list.append(str(ref_seq.id))

                    return '-'.join(new_list)
            elif ordered_list_of_reference_sequences_in_footprint is None:
                # Then this is the final naming of the type which will be done at the end of the analysis
                # once we have given names to all of the reference sequences
                # I have chosen to write this code spearate rather than integrate below
                # as I don't want to break the below
                # code by trying to overwrite orderedListOfRefSeqsInFootprint
                orderedlistofrefseqsinfootprint = self.get_ordered_footprint_list()
                if self.coDom:
                    list_of_majority_reference_sequences = ReferenceSequence.objects.filter(
                        id__in=[int(x) for x in self.MajRefSeqSet.split(',')])
                    # Start the name with the coDom intras in order of abundance.
                    # Then append the noncoDom intras in order of abundance
                    name = '/'.join([
                        str(ref_seq) for ref_seq in orderedlistofrefseqsinfootprint
                        if ref_seq in list_of_majority_reference_sequences])

                    list_of_remaining_reference_sequences = [
                        str(ref_seq) for ref_seq in orderedlistofrefseqsinfootprint
                        if ref_seq not in list_of_majority_reference_sequences]

                    if list_of_remaining_reference_sequences:
                        name += '-{}'.format('-'.join(list_of_remaining_reference_sequences))

                    return name
                else:
                    return '-'.join(str(ref_seq) for ref_seq in orderedlistofrefseqsinfootprint)

            else:
                if self.coDom:
                    list_of_majority_reference_sequences = ReferenceSequence.objects.filter(
                        id__in=[int(x) for x in self.MajRefSeqSet.split(',')])

                    # Start the name with the coDom intras in order of abundance.
                    # Then append the noncoDom intras in order of abundance
                    name = '/'.join(
                        [str(ref_seq) for ref_seq in ordered_list_of_reference_sequences_in_footprint if ref_seq in
                         list_of_majority_reference_sequences])

                    list_of_remaining_reference_sequences = [
                        str(ref_seq) for ref_seq in ordered_list_of_reference_sequences_in_footprint
                        if ref_seq not in list_of_majority_reference_sequences]

                    if list_of_remaining_reference_sequences:
                        name += '-{}'.format('-'.join(list_of_remaining_reference_sequences))
                    return name
                else:
                    return '-'.join(str(ref_seq) for ref_seq in ordered_list_of_reference_sequences_in_footprint)
        return

    def generate_max_min_ratios(self, seqeuence_abundances_for_divs_for_each_clade_collection, ordered_footprint_list):
        # This is the new version using the new concept of lowering the lower limit
        # down to 0.001% if any of the intras are found at below 5% as
        # this indicates that this intra probably spans the withinCladeCutoff
        max_min_list = []
        for j in range(len(ordered_footprint_list)):  # for each refSeqabundance/column
            max_val = 0
            min_val = 1
            # We would use a tuple to hold max min but they are converted to lists upon json.dumps
            max_min_list.append([max_val, min_val])
            # for each row/cladecollection
            for i in range(len(seqeuence_abundances_for_divs_for_each_clade_collection)):
                sample_total_sequences_of_intras_in_type = sum(
                    seqeuence_abundances_for_divs_for_each_clade_collection[i])
                # abundance of intra in q / totalseqs in type for clade_collection_object
                maxmin_candidate = (seqeuence_abundances_for_divs_for_each_clade_collection[i][j]
                                    / sample_total_sequences_of_intras_in_type)

                if maxmin_candidate > max_val:
                    max_val = maxmin_candidate
                    max_min_list[-1][0] = max_val  # Populate the tuple for this refseq with the new max
                if maxmin_candidate < min_val:
                    # If we find any of the intras at a rel abund of < 5% then we lower the lower limit to 0.5% (0.005)
                    # This way we should hopefully mitigate any artefacts caused by the within cladeCutOff
                    # Additionally we need to note that our 0.03 within clade cutoff may be acting on this intra
                    # to cause an artefact. To do this we need to add the intra (ref_seq.id) to the types
                    # artefactIntras list

                    if maxmin_candidate < 0.06:
                        min_val = 0.0001
                        max_min_list[-1][1] = min_val
                        self.note_artefact_intra(ordered_footprint_list[j].id)
                    else:
                        min_val = maxmin_candidate
                        max_min_list[-1][1] = min_val  # Populate the tuple for this refseq with the new min
        self.maxMinRatios = json.dumps(max_min_list)

    def note_artefact_intra(self, reference_sequence_uid):
        # Add the ID of the refseq in question to the artefactIntras list of the type if the ID not already in list
        if self.artefactIntras != '':
            artefact_list = self.artefactIntras.split(',')
            if str(reference_sequence_uid) in artefact_list:
                return
            else:
                artefact_list.append(str(reference_sequence_uid))
                self.artefactIntras = ','.join(artefact_list)
                self.save()
        else:
            self.artefactIntras = str(reference_sequence_uid)
            self.save()

    def get_clade_collections_found_in_initially(self):
        return CladeCollection.objects.filter(
            id__in=[int(x) for x in self.listOfCladeCollectionsFoundInInitially.split(',')])

    def get_clade_collections(self):
        return CladeCollection.objects.filter(id__in=[int(x) for x in self.listOfCladeCollections.split(',')])

    def get_max_min_ratios(self):
        return json.loads(self.maxMinRatios)

    def get_ordered_footprint_list(self):
        # Important that we keep the order of the list
        list_of_ref_seq_uids = [int(a) for a in self.orderedFootprintList.split(',')]
        return [ReferenceSequence.objects.get(id=uid) for uid in list_of_ref_seq_uids]

    def set_maj_ref_seq_set(self, set_of_reference_sequences):
        self.MajRefSeqSet = ','.join(str(refseq.id) for refseq in list(set_of_reference_sequences))

    def get_majority_reference_sequence_set(self):
        list_of_reference_sequence_uids = [int(a) for a in self.MajRefSeqSet.split(',')]
        reference_sequence_query_set = ReferenceSequence.objects.filter(id__in=list_of_reference_sequence_uids)
        return set([refseq for refseq in reference_sequence_query_set])

    def __str__(self):
        return self.name


class CladeCollectionType(models.Model):
    analysisTypeOf = models.ForeignKey(AnalysisType, on_delete=models.CASCADE, null=True)
    # analysisTypeOf = models.IntegerField(null=True)
    cladeCollectionFoundIn = models.ForeignKey(CladeCollection, on_delete=models.CASCADE, null=True)
    # cladeCollectionFoundIn = models.IntegerField(null=True)

    def __str__(self):
        # return ','.join([str(refseq.id) for refseq in self.analysisTypeOf.getOrderedFootprintList()])
        return self.analysisTypeOf.name


# TODO 08/12/17 Consider adding a basal seq of argument in this. This will improve our ability to separate basal
# ITS2 types in type discovery. We would call something basal from a given sequence if it was an exact match
# except for a one bp differentiation. Or I guess we could just go by the closest match. Need to consider I guess what
# would happen to some of the biguns like C21 and C39. But this could be really good.
class ReferenceSequence(models.Model):
    name = models.CharField(max_length=30, default='noName')
    hasName = models.BooleanField(default=False)
    clade = models.CharField(max_length=30)
    sequence = models.CharField(max_length=500)
    accession = models.CharField(max_length=50, null=True)

    def __str__(self):
        if self.hasName:
            return self.name
        else:
            return '{}'.format(self.id)


class DataSetSampleSequence(models.Model):
    cladeCollectionTwoFoundIn = models.ForeignKey(CladeCollection, on_delete=models.CASCADE, null=True)
    referenceSequenceOf = models.ForeignKey(ReferenceSequence, on_delete=models.CASCADE, null=True)
    # referenceSequenceOf = models.IntegerField(null=True)
    abundance = models.IntegerField(default=0)
    data_set_sample_from = models.ForeignKey(DataSetSample, on_delete=models.CASCADE, null=True)

    def __str__(self):
        if self.referenceSequenceOf.hasName:
            return self.referenceSequenceOf.name
        else:
            return 'ID=' + str(self.id)
