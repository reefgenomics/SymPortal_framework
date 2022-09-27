from django.db import models
import json
import general
from datetime import datetime
import sp_config
if sp_config.system_type == 'remote':
    from werkzeug.security import generate_password_hash, check_password_hash

# python3 manage.py graph_models -a -g -o my_project.svg
# You can visualise these models using the following commands in the terminal
# http://cheng.logdown.com/posts/2015/07/07/visualize-the-relationships-of-django-models-using-django-extensions
# sudo apt install graphviz I had to install that in the terminal as well (normal terminal) i.e. not with python3 etc.


class DataSet(models.Model):
    objects = models.Manager()
    name = models.CharField(max_length=60, default='something')
    reference_fasta_database_used = models.CharField(max_length=60, default='None')
    submitting_user = models.CharField(max_length=100, default='no_user_defined')
    submitting_user_email = models.CharField(max_length=100, default='no_email_defined')
    working_directory = models.CharField(max_length=300, default='None')
    time_stamp = models.CharField(max_length=100, default='None')
    loading_complete_time_stamp = models.CharField(max_length=100, default='None')

    def __str__(self):
        return self.name


class DataSetSample(models.Model):

    objects = models.Manager()
    data_submission_from = models.ForeignKey(DataSet, on_delete=models.CASCADE, null=True)
    name = models.CharField(max_length=200, default='None')

    # This is the absolute number of sequences after make.contigs
    num_contigs = models.IntegerField(default=0)
    fastq_fwd_file_name = models.CharField(max_length=100, null=True)
    fastq_fwd_file_hash = models.CharField(max_length=100, null=True)
    fastq_rev_file_name = models.CharField(max_length=100, null=True)
    fastq_rev_file_hash = models.CharField(max_length=100, null=True)

    # store the aboslute number of sequences after inital mothur QC i.e. before tax and size screening
    post_qc_absolute_num_seqs = models.IntegerField(default=0)
    # This is the unique number of sequences after inital mothur QC i.e. before tax and size screening
    post_qc_unique_num_seqs = models.IntegerField(default=0)

    # Absolute number of sequences after sequencing QC and screening for Symbiodiniaceae (i.e. Symbiodiniaceae only)
    absolute_num_sym_seqs = models.IntegerField(default=0)
    # Same as above but the number of unique seqs
    unique_num_sym_seqs = models.IntegerField(default=0)

    # store the abosolute number of sequenes that were not considered Symbiodiniaceae
    non_sym_absolute_num_seqs = models.IntegerField(default=0)
    # This is the number of unique sequences that were not considered Symbiodiniaceae
    non_sym_unique_num_seqs = models.IntegerField(default=0)

    # store the abosulte number of sequences that were lost during the size selection
    size_violation_absolute = models.IntegerField(default=0)
    # store the unique number of sequences that were lost during the size screening
    size_violation_unique = models.IntegerField(default=0)

    # store the number of absolute sequences remaining after MED
    post_med_absolute = models.IntegerField(default=0)
    # store the number of unique sequences remaining after MED (nodes)
    post_med_unique = models.IntegerField(default=0)

    error_in_processing = models.BooleanField(default=False)
    error_reason = models.CharField(max_length=100, default='noError')
    cladal_seq_totals = models.CharField(max_length=5000, null=True)

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

def get_creation_time_stamp_string():
    return str(datetime.now()).split('.')[0].replace('-','').replace(' ','T').replace(':','')

class DataAnalysis(models.Model):
    # This will be a jsoned list of uids of the dataSubmissions that are included in this analysis
    objects = models.Manager()
    list_of_data_set_uids = models.CharField(max_length=5000, null=True)
    within_clade_cutoff = models.FloatField(default=0.04)
    name = models.CharField(max_length=100, null=True)
    description = models.CharField(max_length=5000, null=True)
    time_stamp = models.CharField(max_length=100, default='None')
    submitting_user = models.CharField(max_length=100, default='no_user_defined')
    submitting_user_email = models.CharField(max_length=100, default='no_email_defined')
    analysis_complete_time_stamp = models.CharField(max_length=100, default='None')

    def get_clade_collections(self):
        list_of_uids = [int(x) for x in self.list_of_data_set_uids.split(',')]
        clade_collections = []
        for uid_list in general.chunks(list_of_uids):
            clade_collections.extend(list(CladeCollection.objects.filter(data_set_sample_from__data_submission_from__in=uid_list)))
        return clade_collections

class Study(models.Model):
    objects = models.Manager()
    data_set_samples = models.ManyToManyField(DataSetSample)
    name = models.CharField(max_length=100, null=False, unique=True)
    title = models.CharField(max_length=250, null=True)
    is_published = models.BooleanField(default=False)
    location = models.CharField(max_length=50, null=True)
    run_type = models.CharField(max_length=50, default="remote")
    article_url = models.CharField(max_length=250, null=True)
    data_url = models.CharField(max_length=250, null=True)
    data_explorer = models.BooleanField(default=False)
    display_online = models.BooleanField(default=False)
    analysis = models.BooleanField(default=True)
    data_analysis = models.ForeignKey(DataAnalysis, on_delete=models.SET_NULL, null=True)
    author_list_string = models.CharField(max_length=500, null=True)
    additional_markers = models.CharField(max_length=200, null=True)
    creation_time_stamp = models.CharField(
        max_length=100,
        default=get_creation_time_stamp_string
    )

    def __print__(self):
        print(f'< Study: id {self.id}, name {self.name} >')
        print(f'title: {self.title}')
        print(f'is_published: {self.is_published}')
        print(f'location: {self.location}')
        print(f'run_type: {self.run_type}')
        print(f'article_url: {self.article_url}')
        print(f'data_url: {self.data_url}')
        print(f'data_explorer: {self.data_explorer}')
        print(f'display_online: {self.display_online}')
        print(f'analysis: {self.analysis}')
        print(f'author_list_string: {self.author_list_string}')
        print(f'additional_markers: {self.additional_markers}')
        print(f'creation_time_stamp: {self.creation_time_stamp}')
        print(f'related to {len(list(self.data_set_samples))} '
              f'DataSetSamples objects, the first being {list(self.data_set_samples)[0].name}')

    def __str__(self):
        return f'< Study: id {self.id}, name {self.name} >'

    def __repr__(self):
        return f'< Study: id {self.id}, name {self.name} >'

class Citation(models.Model):
    """
    An object to keep track of SymPortal citations that we are aware of that should not be linked to a Study
    """
    objects = models.Manager()
    title = models.CharField(max_length=250, null=True, unique=True)
    article_url = models.CharField(max_length=250, null=True)
    author_list_string = models.CharField(max_length=500, null=True)
    year = models.CharField(max_length=4, null=True)

class User(models.Model):
    objects = models.Manager()
    name = models.CharField(max_length=100, null=False, unique=True)
    studies = models.ManyToManyField(Study)
    password_hash = models.CharField(max_length=200, null=True)
    is_admin = models.BooleanField(default=False)
    has_upload_permission = models.BooleanField(default=False)

    def __str__(self):
        return f'< User: id {self.id}, name {self.name} >'

    def __repr__(self):
        return f'< User: id {self.id}, name {self.name} >'

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)


class Submission(models.Model):
    """
    A class for representing a user dataset submission.
    It will hold information that will allow the cron jobs to process
    data that has been uploaded to the SymPortal.org webpage
    """
    objects = models.Manager()
    # This name will be used for the DataSet object and the Study object that will be associated to this object
    name = models.CharField(max_length=60, null=False, unique=True)
    # The optional title that will be given to the Study object that is created
    title = models.CharField(max_length=250, null=True)
    # The optional location to be transferred to the associated Study object
    location = models.CharField(max_length=50, null=True)
    # The location of the directory holding the seq files and datasheet on the web hosting server, i.e. linode
    web_local_dir_path = models.CharField(max_length=300, null=False, unique=True)
    # The location of the directory holding the seq files and datasheet on the symportal framework server, i.e. zygote
    framework_local_dir_path = models.CharField(max_length=300, null=False, unique=True)
    # The progress of the submission
    progress_status = models.CharField(max_length=50, null=False, default='pending_submission')
    # If an Error has occured
    error_has_occured = models.BooleanField(default=False)
    # associated DataSet
    associated_dataset = models.ForeignKey(DataSet, on_delete=models.CASCADE, null=True)
    # associated Study
    associated_study = models.ForeignKey(Study, on_delete=models.CASCADE, null=True)
    # submitting User
    submitting_user = models.ForeignKey(User, on_delete=models.CASCADE, null=False)
    # number of samples
    number_samples = models.IntegerField(null=False, default=0)
    # These fields will hold the times that various checkpoints are reached
    # submission time
    submission_date_time = models.CharField(max_length=25, null=False, default=get_creation_time_stamp_string)
    transfer_to_framework_server_date_time = models.CharField(max_length=25, null=True)
    loading_started_date_time = models.CharField(max_length=25, null=True)
    loading_complete_date_time = models.CharField(max_length=25, null=True)
    analysis_started_date_time = models.CharField(max_length=25, null=True)
    analysis_complete_date_time = models.CharField(max_length=25, null=True)
    study_output_started_date_time = models.CharField(max_length=25, null=True)
    study_output_complete_date_time = models.CharField(max_length=25, null=True)
    transfer_to_web_server_date_time = models.CharField(max_length=25, null=True)
    # Whether this submission should be going into an analysis or not
    # I.e. if it contains seawater samples or something similar then it should not be going into analysis
    # conservatively set as False
    for_analysis = models.BooleanField(default=False, null=False)

    # The path to the directory in which the result files are output for:
    # Framework server
    framework_results_dir_path = models.CharField(max_length=300, null=True)
    # Web server
    web_results_dir_path = models.CharField(max_length=300, null=True)

class CladeCollection(models.Model):
    objects = models.Manager()
    data_set_sample_from = models.ForeignKey(DataSetSample, on_delete=models.CASCADE, null=True)
    clade = models.CharField(max_length=1)
    # the method below to get the footprint of the clade_collection_object is incredibly slow.
    # Hopefully a much faster way will be to store the ID's of the refseqs that make up the footprint
    # I will therefore store the reference uids in the field below
    footprint = models.CharField(max_length=100000, default=True)

    # This will return the foot print of the analysedSampleSequences that are found above the given percentage cutoff
    def cutoff_footprint(self, cutoff):
        # get total seqs in cladeCollection
        total = 0

        for dsss in DataSetSampleSequence.objects.filter(clade_collection_found_in=self):
            total += dsss.abundance
        sequence_number_cutoff = cutoff * total
        frset = set([])
        for dsss in DataSetSampleSequence.objects.filter(clade_collection_found_in=self):
            if dsss.abundance > sequence_number_cutoff:
                frset.add(dsss.reference_sequence_of)
        return frozenset(frset)

    def __str__(self):
        return self.data_set_sample_from.name


class AnalysisType(models.Model):
    objects = models.Manager()
    data_analysis_from = models.ForeignKey(DataAnalysis, on_delete=models.CASCADE, null=True)
    # This should be a frozen set of referenceSequences
    # As this is not possible in a django field or jsonable
    # I will instead make it a string of reference sequence uids
    # I will make set and get methods for the footprint
    # This will be a list of refSeqs that make up the footprint in order of their abundance when type first defined
    # This is a commar separated string of the uids of the ref seqs that define the type
    ordered_footprint_list = models.CharField(max_length=200, null=True)
    # Same for listOfMajs
    # set() of refseqs that are Majs in each of the CCs this type was initially identified in.
    # Note that this is therefore in no particular order
    majority_reference_sequence_set = models.CharField(max_length=40, null=True)

    list_of_clade_collections = models.CharField(max_length=1000000, null=True)
    # This is a 2D list, a list for each clade collection in order of the listofCladeCollections
    # Within each list the absolute abundances of the defining seqs in order of ordered_footprint_list
    footprint_sequence_abundances = models.CharField(max_length=1000000, null=True)

    # Same as above but the proportion of the seqs to each other in the cladecollection.
    footprint_sequence_ratios = models.CharField(max_length=1000000, null=True)

    clade = models.CharField(max_length=1)
    co_dominant = models.BooleanField(default=False)

    name = models.CharField(max_length=1000, null=True)

    # The list of speceis that this type is associated with
    species = models.CharField(max_length=200, null=True)

    # this list will keep track of which of the defining intras of this type are 'unlocked' i.e. at least one
    # of the instances of that intra were found at <5%. We will use this list to the 'artefact type creation'
    # This artefact_intras will hold a char string of comma separated ints that represent the id's of the
    # refseqs that are unlocked
    artefact_intras = models.CharField(max_length=5000, default='')

    def get_ratio_list(self):
        return json.loads(self.footprint_sequence_ratios)

    def get_clade_collections(self):
        uids_for_query = [int(x) for x in self.list_of_clade_collections.split(',')]
        cc_objs_list = []
        for uid_list in general.chunks(uids_for_query):
            cc_objs_list.extend(list(CladeCollection.objects.filter(id__in=uid_list)))
        return cc_objs_list

    def __str__(self):
        return self.name


class CladeCollectionType(models.Model):
    objects = models.Manager()
    analysis_type_of = models.ForeignKey(AnalysisType, on_delete=models.CASCADE, null=True)
    # analysis_type_of = models.IntegerField(null=True)
    clade_collection_found_in = models.ForeignKey(CladeCollection, on_delete=models.CASCADE, null=True)
    # clade_collection_found_in = models.IntegerField(null=True)

    def __str__(self):
        # return ','.join([str(refseq.id) for refseq in self.analysis_type_of.get_ordered_footprint_list()])
        return self.analysis_type_of.name


class ReferenceSequence(models.Model):
    objects = models.Manager()
    name = models.CharField(max_length=30, default='noName')
    has_name = models.BooleanField(default=False)
    clade = models.CharField(max_length=30)
    sequence = models.CharField(max_length=500)
    accession = models.CharField(max_length=50, null=True)

    def __str__(self):
        if self.has_name:
            return self.name
        else:
            return f'{self.id}_{self.clade}'


class DataSetSampleSequence(models.Model):
    objects = models.Manager()
    clade_collection_found_in = models.ForeignKey(CladeCollection, on_delete=models.CASCADE, null=True)
    reference_sequence_of = models.ForeignKey(ReferenceSequence, on_delete=models.CASCADE, null=True)
    # reference_sequence_of = models.IntegerField(null=True)
    abundance = models.IntegerField(default=0)
    data_set_sample_from = models.ForeignKey(DataSetSample, on_delete=models.CASCADE, null=True)

    def __str__(self):
        if self.reference_sequence_of.has_name:
            return self.reference_sequence_of.name
        else:
            return 'ID=' + str(self.id)

class DataSetSampleSequencePM(models.Model):
    # this is the pre-MED version of the DataSetSampleSequence object
    # its purpose is to keep track of the
    objects = models.Manager()
    data_set_sample_from = models.ForeignKey(DataSetSample, on_delete=models.CASCADE, null=True)
    reference_sequence_of = models.ForeignKey(ReferenceSequence, on_delete=models.CASCADE, null=True)
    # reference_sequence_of = models.IntegerField(null=True)
    abundance = models.IntegerField(default=0)


    def __str__(self):
        if self.reference_sequence_of.has_name:
            return self.reference_sequence_of.name
        else:
            return 'ID=' + str(self.id)
