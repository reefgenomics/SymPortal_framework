from django.db import models
from django.db.models import Max
from datetime import datetime
import os
import json

#python3 manage.py graph_models -a -g -o my_project.svg
# You can visualise these models using the following commands in the terminal
#http://cheng.logdown.com/posts/2015/07/07/visualize-the-relationships-of-django-models-using-django-extensions
## sudo apt install graphviz I had to install that in the terminal as well (normal terminal) i.e. not with python3 etc.


class symportal_framework(models.Model):
    latest_reference_fasta = models.CharField(max_length = 30, default='symClade_2_2.fa')
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

class data_set(models.Model):
    name = models.CharField(max_length = 60, default='something')
    reference_fasta_database_used = models.CharField(max_length = 60, default='None')
    submittingUser = models.CharField(max_length=30, default='Bob')
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
    #http://stackoverflow.com/questions/28896627/executing-code-on-startup-in-django-1-7
    currentlyBeingProcessed = models.BooleanField(default=False)
    timeStamp = models.CharField(max_length=40, default='None')
    #http://stackoverflow.com/questions/27220480/django-datetime-migration-error
    # We can programatically set whenSubmitted by '= datetime.now()'
    # The Datetimes were causing a lot of problems when doing the dumpdata and loaddata. I think it was becuase of the
    # for mat of the date times, specifically because they had a time zone applied.
    # I think we can likely fix this by specifying formats for the datetimes in the settings file
    # We can then make sure thatn when we enter a datetime using the python that we use the correct format.
    # https://docs.djangoproject.com/en/1.9/ref/settings/#datetime-input-formats


    def __str__(self):
        return self.name

    #http://stackoverflow.com/questions/11210396/django-set-field-values-dynamically-based-on-other-fields
    #http://stackoverflow.com/questions/4380879/django-model-field-default-based-off-another-field-in-same-model
    def save(self, *args, **kwargs):
        # This used to have content for percent complete to be calculated but I have gotten rid of this
        # I have left this structure here in case we want to use it again.
        # Content would go here
        super().save(*args, **kwargs)


class data_set_sample(models.Model):

    dataSubmissionFrom = models.ForeignKey(data_set, on_delete=models.CASCADE, null=True)
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
    #This is the number of unique sequences that were not considered Symbiodinium
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

    def __str__(self):
        return self.name


class data_analysis(models.Model):
    # This will be a jsoned list of IDs of the dataSubmissions that are included in this analysis
    listOfDataSubmissions = models.CharField(max_length=100, null=True)
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
    timeStamp = models.CharField(max_length=40, default='None')
    def getCladeCollections(self):
        listOfIds = [int(x) for x in self.listOfDataSubmissions.split(',')]
        cladeCollections = clade_collection.objects.filter(dataSetSampleFrom__dataSubmissionFrom__in=listOfIds)
        return cladeCollections


class clade_collection(models.Model):
    dataSetSampleFrom = models.ForeignKey(data_set_sample, on_delete=models.CASCADE, null=True)
    clade = models.CharField(max_length=1)
    #the method below to get the footprint of the CC is incredibly slow.
    # Hopefully a much faster way will be to store the ID's of the refseqs that make up the footprint
    # I will therefore store the reference IDs in the field below
    footPrint = models.CharField(max_length=100000, default=True)



    # maj() should return the analysedSampleSequence associated with the
    # cladeCollection in Q with the largest abundance
    def maj(self):
        maxAbundance = \
        data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=self).aggregate(Max('abundance'))[
            'abundance__max']
        # A .get() could fail if there are two sequences with an equal abundance so I
        # will try the get first and then failing that I will do a .filter() and return the first element
        try:
            return data_set_sample_sequence.objects.get(cladeCollectionTwoFoundIn=self, abundance=maxAbundance)

        except:
            return data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=self, abundance=maxAbundance)[0]

    # range will be how many of the seqs we want to return
    def orderedListOfRefSeqs(self, all):
        if all:
            dsss = list(data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=self).order_by('-abundance'))
            orderedRefSeqs = [ds.referenceSequenceOf for ds in dsss]
        else:
            # Try is for incase there are not 10 seqs in the CC
            try:
                dsss = list(data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=self).order_by('-abundance')[:10])
                orderedRefSeqs = [ds.referenceSequenceOf for ds in dsss]
            except:
                dsss = list(data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=self).order_by('-abundance'))
                orderedRefSeqs = [ds.referenceSequenceOf for ds in dsss]
        return orderedRefSeqs

    # This will return the foot print of the analysedSampleSequences that are found above the given percentage cutoff
    def cutOffFootprint(self, cutoff):
        # get total seqs in cladeCollection
        total = 0
        dsss = list(data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=self))

        for dsss in data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=self):
            total += dsss.abundance
        seqNumCutoff = cutoff * total
        frset = set([])
        for dsss in data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=self):
            if dsss.abundance > seqNumCutoff:
                frset.add(dsss.referenceSequenceOf)
        return frozenset(frset)

    # footprint should return a frozen of referenceSequenceObjects associated
    # with the analysedSampleSequences associated with the cladeCollectionInQ
    def footprint(self):
        return frozenset(refseq for refseq in reference_sequence.objects.filter(id__in=[int(a) for a in self.footPrint.split(',')]))

    def __str__(self):
        return self.dataSetSampleFrom.name


class analysis_group(models.Model):
    dataAnalysisFrom = models.ForeignKey(data_analysis, on_delete=models.CASCADE, null=True)
    name = models.CharField(max_length=100, null=True)

    # We will generate a name for the analysis_group at some point later according to the analysisTypes it contains


class analysis_type(models.Model):
    dataAnalysisFrom = models.ForeignKey(data_analysis, on_delete=models.CASCADE, null=True)
    analysisGroupOf = models.ForeignKey(analysis_group, on_delete=models.SET_NULL, null=True)
    # This should be a frozen set of referenceSequences
    # As this is not possible in a django field or jsonable
    # I will instead make it a string of reference sequence IDs
    # I will make set and get methods for the footprint
    # This will be a list of refSeqs that make up the footprint in order of their abundance when type first defined
    # This is a commar separated string of the IDs of the ref seqs that define the type
    orderedFootprintList = models.CharField(max_length=200, null=True)
    # Same for listOfMajs
    # set() of refseqs that are Majs in each of the CCs this type was initially identified in.
    # Note that this is therefore in no particular order
    MajRefSeqSet = models.CharField(max_length=40, null=True)
    # Same for this
    # The list of cladeCollections in which the type was defined from
    listOfCladeCollectionsFoundInInitially = models.CharField(max_length=5000, null=True)
    # The list of cladeCollections that the type was associated with after one iteration of assigningTypes
    listOfCladeCollections = models.CharField(max_length=10000, null=True)
    # This is a 2D list, a list for each clade collection in order of the listofCladeCollections
    # Within each list the absolute abundances of the defining seqs in order of orderedFootprintList
    footprintSeqAbundances = models.CharField(max_length=10000, null=True)
    # Same as above but the relative abundance of the seq in Q as a function of all of the sequences in the cladeCollection
    footprintSeqRatios = models.CharField(max_length=100000, null=True)
    clade = models.CharField(max_length=1)
    coDom = models.BooleanField(default=False)

    name = models.CharField(max_length=1000, null=True)
    # [(max,min) for refseq in ordered footprint list]
    maxMinRatios = models.CharField(max_length=10000, null=True)
    # The list of speceis that this type is associated with
    species = models.CharField(max_length=200, null=True)

    #this list will keep track of which of the defining intras of this type are 'unlocked' i.e. at least one
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

    def initTypeAttributes(self, listOfCC, footprintlistofrefseqs):
        # Becuase this can be run twice on a type due to the checking for within clade cutoff methods
        # We need to reset the artefactIntras
        self.artefactIntras = ''

        ### DENBUG ###
        if list(footprintlistofrefseqs)[0].name == 'C3cf':
            foo = 'bar'


        #listOfCC = [seq.cladeCollectionFoundIn for seq in listOfMajSeqs]

        try:
            self.listOfCladeCollectionsFoundInInitially = ','.join([str(cc.id) for cc in list(listOfCC)])
        except:
            this = footprintlistofrefseqs
            apples = 'fo'
        footprintlist = list(footprintlistofrefseqs)
        countList = [[] for cc in listOfCC]

        # Here get the abundance of each of the defining intras in each of the samples of this type
        # Then we can sort it and make a name

        # TODO here we are going to introduce a new concept.
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
        # If at least one of the intras in one of the CC in the type is found, then we will consider this profile
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
        setOfMajRefSeqs = set()
        for CC in listOfCC:
            maxAbund = 0
            majRefSeq = None
            listOfSeqsFromCC = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=CC)
            for refSeq in footprintlistofrefseqs:
                sampleSeqOfRefSeqType = listOfSeqsFromCC.get(referenceSequenceOf=refSeq)
                abundanceOfRefSeqInQ = sampleSeqOfRefSeqType.abundance
                if abundanceOfRefSeqInQ > maxAbund:
                    majRefSeq = sampleSeqOfRefSeqType.referenceSequenceOf
                    maxAbund = sampleSeqOfRefSeqType.abundance
                countList[listOfCC.index(CC)].append(abundanceOfRefSeqInQ)
            setOfMajRefSeqs.add(majRefSeq)

        # if setOfMajRefSeqs > 1 then this is coDom
        if len(setOfMajRefSeqs) > 1:
            self.coDom = True
        else:
            self.coDom = False
        self.setMajRefSeqSet(setOfMajRefSeqs)
        # At this point we have the abundance of each sequence for each instance of the type
        # we have 2d array where list are cladeCollections, and items in each list represent sequence
        # Now create counting list to count the abundance of each seq in the footprints refseqs
        try:
            abundanceList = [[refseq, 0] for refseq in footprintlist]
            for i in range(len(countList)): # for each cladeCollection with footprint
                for j in range(len(countList[0])): # for each refseq in footprint
                    abundanceList[j][1] += countList[i][j]
        except:
            foo = 'bar'
        #order the list
        orderedFootprintList = [a[0] for a in sorted(abundanceList, key=lambda x: x[1], reverse=True)]
        # This now gives us the abundance order of the refseqs
        # TODO 08/12/17 this is currently working with absolute abundances instead of relative abundances
        # Which strikes me as incorrect. BUT I don't want to start messing with this now.
        # just so long as we work with the relative abundances as funtion of the type seuences in the CC


        self.orderedFootprintList = ','.join([str(refseq.id) for refseq in orderedFootprintList])

        if not self.isLockedType:
            self.name = self.generateName(orderedFootprintList)

        # The below code re orders the countlist 2D list so that the columns are in order of refseq abundance
        countListOrdered = []
        for row in countList:
            countListOrdered.append([]) # new empty list for each row
            for element in orderedFootprintList: #for each refseq in order of abundance
                countListOrdered[-1].append(row[footprintlist.index(element)])

        # recover the memory from countList
        del countList
        # Here we have the 2d list of raw abundances of seqs for every clade collection (row) for every seq (column)
        # Now we can store it in a char field by json.dumpsing it.
        self.footprintSeqAbundances = json.dumps(countListOrdered)

        # These now represent relative abundances rather than ratios

        self.generateMaxMinRatios(countListOrdered, orderedFootprintList)
        self.generateRatioList()
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

    def updateTypeAttributes(self):
        # Important to use this get() format rather than all().filter() format as we need to maintain the order of the cladeCollection list
        # NB this is now updating from listoflcadecollections rather than from listoflcadecollectionsfoundininitially

        #### DEBUG ####
        if self.name == 'A4-A4b-A4c':
            pausey = 'pausey'

        try:
            listOfCC = [clade_collection.objects.get(id=ID) for ID in [int(x) for x in self.listOfCladeCollections.split(',')]]
        except:
            apples = 'asdf'
        orderedListOfRefSeqsInType = [reference_sequence.objects.get(id=ID) for ID in [int(x) for x in self.orderedFootprintList.split(',')]]

        countList = [[] for cc in listOfCC]

        for CC in listOfCC:
            listOfSeqsFromCC = data_set_sample_sequence.objects.filter(cladeCollectionTwoFoundIn=CC)
            for refSeq in orderedListOfRefSeqsInType:
                sampleSeqOfRefSeqType = listOfSeqsFromCC.get(referenceSequenceOf=refSeq)
                countList[listOfCC.index(CC)].append(sampleSeqOfRefSeqType.abundance)

        #Here I will update/initiate self.MajRefSeqSet and self.coDom using the countList info
        MajRefSeqSet = set()
        for CClist in countList:
            try:
                maxValue = max(CClist)
            except:
                pasd = 'asdf'
            index = CClist.index(maxValue)
            MajRefSeqSet.add(orderedListOfRefSeqsInType[index])
        self.MajRefSeqSet = ','.join([str(a.id) for a in list(MajRefSeqSet)])
        if len(MajRefSeqSet) > 1:
            self.coDom = True
        else:
            self.coDom = False
        self.save()
        # At this point we have the abundance of each sequence for each instance of the type
        # we have 2d array where list are samples and items in list represent sequence
        # Now create counting list to count the abundance of each seq in the footprints refseqs
        abundanceList = [[refseq, 0] for refseq in orderedListOfRefSeqsInType]
        for i in range(len(countList)):  # for each cladeCollection with footprint
            for j in range(len(countList[0])):  # for each refseq in footprint
                try:
                    abundanceList[j][1] += countList[i][j]
                except:
                    apples = 'adsf'
        orderedFootprintList = [a[0] for a in sorted(abundanceList, key=lambda x: x[1], reverse=True)]
        # Because this is also used when creating a brand new type in the bimodal distribution splitting
        # We should always create a new name
        # if self.orderedFootprintList != ','.join([str(refseq.id) for refseq in orderedFootprintList]):
            #Then the order of the footprint has changed and the name needs to be redone
            #Else it doesn't need to be redone
        self.orderedFootprintList = ','.join([str(refseq.id) for refseq in orderedFootprintList])
        if self.orderedFootprintList == '':
            apples = 'adsf'

        if not self.isLockedType:
            self.name = self.generateName(orderedFootprintList)
        # If the order of the footPrint hasn't changed then the count list is already in the right order and we don't need to do the code below.
        # If the order has changed then we do need to redo the code
        # The below code re orders the countlist 2D list os that the columns are in order of refseq abundance
        countListOrdered = []
        for row in countList:
            countListOrdered.append([])  # new empty list for each row
            for element in orderedFootprintList:  # for each refseq in order of abundance
                countListOrdered[-1].append(row[orderedListOfRefSeqsInType.index(element)])
        # recover the memory from countList
        del countList

        # Here we have the 2d list of raw abundances of seqs for every clade collection (row) for every seq (column)
        # Now we can store it in a char field by json.dumpsing it.

        self.footprintSeqAbundances = json.dumps(countListOrdered)
        # These now represent relative abundances rather than ratios
        self.generateMaxMinRatios(countListOrdered, orderedFootprintList)
        self.generateRatioList()

    def addCCToCladeCollectionList(self, CC):
        if self.listOfCladeCollections:# If empty then no need to add commar
            self.listOfCladeCollections += ',{}'.format(CC.id)
        else:
            self.listOfCladeCollections = str(CC.id)

    def addCCListToCladeCollectionList(self, CCidList):
        if self.listOfCladeCollections: # If empty then no need to add commar
            self.listOfCladeCollections += ',{}'.format(','.join([str(a) for a in CCidList]))
        else:
            self.listOfCladeCollections = ','.join([str(a) for a in CCidList])

    def removeCCListFromInitialCladeCollectionList(self, CCidList):
        listOfIDsAsStr = self.listOfCladeCollectionsFoundInInitially.split(',')
        newList = [id for id in listOfIDsAsStr if id not in CCidList]
        self.listOfCladeCollectionsFoundInInitially = ','.join(newList)
        self.save()


    def removeCCFromCladeCollectionList(self, CC):
        listOfIDsAsStr = self.listOfCladeCollections.split(',')
        listOfIDsAsStr.remove(str(CC.id))
        self.listOfCladeCollections =  ','.join(listOfIDsAsStr)
        self.save()


    # This is the relative abundance of each sequence in the type as a function of the total sequence found
    # in the clade Collection THAT ARE ALSO in the type. SO proportion of type sequences in cladecollection

    def generateRatioList(self):
        newRatioList = []
        footprintSeqAbundancesLoaded = json.loads(self.footprintSeqAbundances)
        for i in range(len(footprintSeqAbundancesLoaded)):# For each CC the type is associated with
            newRatioList.append([])
            denominator = sum(footprintSeqAbundancesLoaded[i])
            for j in range(len(footprintSeqAbundancesLoaded[i])):# For each intra in the type
                ratio = float("%.3f" % (footprintSeqAbundancesLoaded[i][j]/denominator))
                newRatioList[-1].append(ratio)
        self.footprintSeqRatios = json.dumps(newRatioList)

    def getRatioList(self):
        return json.loads(self.footprintSeqRatios)

    def generateName(self, orderedListOfRefSeqsInFootprint=None, accession=None):
        # If self.isLockedType then we skip the naming.
        if not self.isLockedType:

            #this is interesting, sometimes one of the intras that is not a coDom is more abundant than one
            # of the coDom intras to solve this problem we should add the coDoms intas to the name first and then
            # go in the remaining order of the orderedListOfRefSeqsInFootprint

            if accession:
                # Here we output the accession name for the type in Q
                # This checks to see if there are accession numbers available for each of the refseqs in the type
                # If there are accessions available we use these, in the same / and - format as the usual name
                # Where there aren't accessions available for the ref seqs we use the ref seq's ID.
                orderedlistofrefseqsinfootprint = self.getOrderedFootprintList()
                if self.coDom:
                    listOfMajRefSeqs = reference_sequence.objects.filter(
                        id__in=[int(x) for x in self.MajRefSeqSet.split(',')])
                    # Start the name with the coDom intras in order of abundance.
                    # Then append the noncoDom intras in order of abundance
                    orderedListOfMajRefSeqs = [refSeq for refSeq in orderedlistofrefseqsinfootprint if
                                               refSeq in listOfMajRefSeqs]

                    # temp_list = []
                    # for refSeq in orderedListOfMajRefSeqs:
                    #     if refSeq.accession != None:
                    #         temp_list.append(refSeq.accession)
                    #     else:
                    #         temp_list.append(str(refSeq.id))
                    # name = '/'.join(temp_list)
                    name = '/'.join(
                        [refSeq.accession if refSeq.accession != None else str(refSeq.id) for refSeq in orderedListOfMajRefSeqs])
                    listOfRemainingRefSeqs = [refSeq for refSeq in orderedlistofrefseqsinfootprint if
                                              refSeq not in listOfMajRefSeqs]
                    if listOfRemainingRefSeqs:
                        name += '-{}'.format('-'.join(
                            [refSeq.accession if refSeq.accession != None else str(refSeq.id) for refSeq in
                             listOfRemainingRefSeqs]))
                    return name
                else:
                    new_list = []
                    for refSeq in orderedlistofrefseqsinfootprint:
                        if refSeq.accession and refSeq.accession != None:
                            new_list.append(refSeq.accession)
                        else:
                            new_list.append(str(refSeq.id))


                    return '-'.join(new_list)
            elif orderedListOfRefSeqsInFootprint == None:
                # Then this is the final naming of the type which will be done at the end of the analysis
                # once we have given names to all of the reference sequences
                # I have chosen to write this code spearate rather than integrate below as I don't want to break the below
                # code by trying to overwrite orderedListOfRefSeqsInFootprint
                orderedlistofrefseqsinfootprint = self.getOrderedFootprintList()
                if self.coDom:
                    listOfMajRefSeqs = reference_sequence.objects.filter(id__in=[int(x) for x in self.MajRefSeqSet.split(',')])
                    # Start the name with the coDom intras in order of abundance.
                    # Then append the noncoDom intras in order of abundance
                    name = '/'.join([str(refSeq) for refSeq in orderedlistofrefseqsinfootprint if refSeq in listOfMajRefSeqs])
                    listOfRemainingRefSeqs = [str(refSeq) for refSeq in orderedlistofrefseqsinfootprint if refSeq not in listOfMajRefSeqs]
                    if listOfRemainingRefSeqs:
                        name += '-{}'.format('-'.join(listOfRemainingRefSeqs))
                    return name
                else:
                    return '-'.join(str(refSeq) for refSeq in orderedlistofrefseqsinfootprint)

            else:
                if self.coDom:
                    listOfMajRefSeqs = reference_sequence.objects.filter(id__in=[int(x) for x in self.MajRefSeqSet.split(',')])

                    #Start the name with the coDom intras in order of abundance.
                    #Then append the noncoDom intras in order of abundance
                    name = '/'.join([str(refSeq) for refSeq in orderedListOfRefSeqsInFootprint if refSeq in listOfMajRefSeqs])
                    listOfRemainingRefSeqs = [str(refSeq) for refSeq in orderedListOfRefSeqsInFootprint if refSeq not in listOfMajRefSeqs]
                    if listOfRemainingRefSeqs:
                        name += '-{}'.format('-'.join(listOfRemainingRefSeqs))
                    return name
                else:
                    return '-'.join(str(refSeq) for refSeq in orderedListOfRefSeqsInFootprint)
        return

    def generateMaxMinRatios(self, seqCountForDefiningIntrasForEachCC, orderedFootprintList):
        # This is the new version using the new concept of lowering the lower limit
        # down to 0.001% if any of the intras are found at below 5% as this indicates that this intra probably spans the withinCladeCutoff
        maxMinList = []
        for j in range(len(orderedFootprintList)):  # for each refSeqabundance/column
            max = 0
            min = 1
            # We would use a tuple to hold max min but they are converted to lists upon json.dumps
            maxMinList.append([max, min])
            for i in range(len(seqCountForDefiningIntrasForEachCC)): # for each row/cladecollection
                sampleTotalSequencesOfIntrasInType = sum(seqCountForDefiningIntrasForEachCC[i])
                maxminCandidate = seqCountForDefiningIntrasForEachCC[i][j]/sampleTotalSequencesOfIntrasInType # abundance of intra in q / totalseqs in type for CC
                if maxminCandidate > max:
                    max = maxminCandidate
                    maxMinList[-1][0] = max # Populate the tuple for this refseq with the new max
                if maxminCandidate < min:
                    # If we find any of the intras at a rel abund of < 5% then we lower the lower limit to 0.5% (0.005)
                    # This way we should hopefully mitigate any artefacts caused by the within cladeCutOff
                    # Additionally we need to note that our 0.03 within clade cutoff may be acting on this intra
                    # to cause an artefact. To do this we need to add the intra (refSeq.id) to the types
                    # artefactIntras list

                    if maxminCandidate < 0.06:
                        min = 0.0001
                        maxMinList[-1][1] = min
                        self.noteArtefactIntra(orderedFootprintList[j].id)
                    else:
                        min = maxminCandidate
                        maxMinList[-1][1] = min  # Populate the tuple for this refseq with the new min
        self.maxMinRatios = json.dumps(maxMinList)

    def noteArtefactIntra(self, refSeqID):
        # Add the ID of the refseq in question to the artefactIntras list of the type if the ID not already in list
        if self.artefactIntras != '':
            artefactList = self.artefactIntras.split(',')
            if str(refSeqID) in artefactList:
                return
            else:
                artefactList.append(str(refSeqID))
                self.artefactIntras = ','.join(artefactList)
                self.save()
        else:
            self.artefactIntras = str(refSeqID)
            self.save()

    def generateMaxMinRatiosOld(self, seqCountForDefiningIntrasForEachCC, orderedFootprintList):
        # This is the original implementation of this. The version above uses the new concept of lowering the lower limit
        # down to 0.05% if any of the intras are found at below 5% as this indicates that this intra probably spans the withinCladeCutoff
        maxMinList = []
        for j in range(len(orderedFootprintList)):  # for each refSeqabundance/column
            max = 0
            min = 1
            # We would use a tuple to hold max min but they are converted to lists upon json.dumps
            maxMinList.append([max, min])
            for i in range(len(seqCountForDefiningIntrasForEachCC)): # for each row/cladecollection
                sampleTotalSequencesOfIntrasInType = sum(seqCountForDefiningIntrasForEachCC[i])
                maxminCandidate = seqCountForDefiningIntrasForEachCC[i][j]/sampleTotalSequencesOfIntrasInType # abundance of intra in q / abundance of maj intra
                if maxminCandidate > max:
                    max = maxminCandidate
                    maxMinList[-1][0] = max # Populate the tuple for this refseq with the new max
                if maxminCandidate < min:
                    min = maxminCandidate
                    maxMinList[-1][1] = min  # Populate the tuple for this refseq with the new min
        self.maxMinRatios = json.dumps(maxMinList)

    def getCladeCollectionsFoundInInitially(self):
        return clade_collection.objects.filter(id__in=[int(x) for x in self.listOfCladeCollectionsFoundInInitially.split(',')])

    def getCladeCollections(self):
        try:
            return clade_collection.objects.filter(id__in=[int(x) for x in self.listOfCladeCollections.split(',')])
        except:
            applie = 'asdf'
        return


    def getMaxMinRatios(self):
        try:
            return json.loads(self.maxMinRatios)
        except:
            apples = 'pears'


    def getOrderedFootprintList(self):
        # Important that we keep the order of the list
        listOfRefSeqIds = [int(a) for a in self.orderedFootprintList.split(',')]
        return [reference_sequence.objects.get(id=id) for id in listOfRefSeqIds]


    def setMajRefSeqSet(self, setOfRefSeqs):
        self.MajRefSeqSet = ','.join(str(refseq.id) for refseq in list(setOfRefSeqs))

    def getMajRefSeqSet(self):
        listOfRefSeqIds = [int(a) for a in self.MajRefSeqSet.split(',')]
        refSeqQuerySet = reference_sequence.objects.filter(id__in=listOfRefSeqIds)
        return set([refseq for refseq in refSeqQuerySet])

    def __str__(self):
        return self.name


class clade_collection_type(models.Model):
    analysisTypeOf = models.ForeignKey(analysis_type, on_delete=models.CASCADE, null=True)
    # analysisTypeOf = models.IntegerField(null=True)
    cladeCollectionFoundIn = models.ForeignKey(clade_collection, on_delete=models.CASCADE, null=True)
    # cladeCollectionFoundIn = models.IntegerField(null=True)

    def __str__(self):
        #return ','.join([str(refseq.id) for refseq in self.analysisTypeOf.getOrderedFootprintList()])
        return self.analysisTypeOf.name


#TODO 08/12/17 Consider adding a basal seq of argument in this. This will improve our ability to separate basal
# ITS2 types in type discovery. We would call something basal from a given sequence if it was an exact match
# except for a one bp differentiation. Or I guess we could just go by the closest match. Need to consider I guess what
# would happen to some of the biguns like C21 and C39. But this could be really good.
class reference_sequence(models.Model):
    name = models.CharField(max_length=30, default='noName')
    hasName = models.BooleanField(default=False)
    clade = models.CharField(max_length=30)
    sequence = models.CharField(max_length=500)
    accession = models.CharField(max_length=50, null=True)
    def __str__(self):
        if self.name == 'noName':
            return '{}'.format(self.id)
        else:
            return self.name


class data_set_sample_sequence(models.Model):
    cladeCollectionTwoFoundIn = models.ForeignKey(clade_collection, on_delete=models.CASCADE, null=True)
    referenceSequenceOf = models.ForeignKey(reference_sequence, on_delete=models.CASCADE, null=True)
    # referenceSequenceOf = models.IntegerField(null=True)
    abundance = models.IntegerField(default=0)
    data_set_sample_from = models.ForeignKey(data_set_sample, on_delete=models.CASCADE, null=True)


    def __str__(self):
        try:
            return self.referenceSequenceOf.name
        except:
            return 'ID=' + self.id






