"""A method that relies on Django feature. E.g. model classes."""
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
from dbApp.models import DataSet, DataAnalysis, DataSetSample, Study, User
import pandas as pd
import sys
from collections import Counter
from numpy import NaN
from django.core.exceptions import ObjectDoesNotExist
from datetime import datetime
from general import check_lat_lon


def delete_data_set(uid):
    DataSet.objects.get(id=uid).delete()


def delete_data_analysis(uid):
    DataAnalysis.objects.get(id=uid).delete()

def write_ref_seq_objects_to_fasta(path, list_of_ref_seq_objs, identifier='name'):
    with open(path, 'w') as f:
        for ref_seq_obj in list_of_ref_seq_objs:
            if identifier == 'name':
                f.write(f'>{ref_seq_obj.name}\n')
            elif identifier == 'id':
                f.write(f'>{ref_seq_obj.id}\n')
            f.write(f'{ref_seq_obj.sequence}\n')

class ApplyDatasheetToDataSetSamples:
    """Class responsible for allowing us to apply a datasheet to a given set of DataSetSample objects
    that belong to a given DataSet. For the time being we will just be working with single DataSets. In the
    future it may be that we can work with sets of DataSet objects or even custom DataSetSample UID sets"""

    def __init__(self, data_set_uid, data_sheet_path):
        self.ds_obj = DataSet.objects.get(id=int(data_set_uid))
        self.datasheet_path = data_sheet_path
        self.sample_meta_info_df = None
        self._make_data_sheet_df()
        self.dss_objects = [dss for dss in DataSetSample.objects.filter(data_submission_from=self.ds_obj) if dss.name in self.sample_meta_info_df.index.values.tolist()]
        self.dss_object_names_list = [dss.name for dss in self.dss_objects]
        self._check_all_samples_found_exit_if_not()
        self.meta_category_titles = ['sample_type', 'host_phylum', 'host_class', 'host_order', 'host_family',
                                     'host_genus', 'host_species', 'collection_latitude', 'collection_longitude',
                                     'collection_date', 'collection_depth']

    def _check_all_samples_found_exit_if_not(self):
        if not len(self.dss_objects) == len(self.sample_meta_info_df.index):
            # If we are missing some samples report on which samples were missing
            for dss_ind_name in self.sample_meta_info_df.index:
                if dss_ind_name not in self.dss_object_names_list:
                    print(
                        f'{dss_ind_name} from your datasheet was not found in dataset {self.ds_obj.name}; id {self.ds_obj.id}.')
            raise RuntimeError('ABORTING application of datasheet to samples. No database objects have been changed')

    def apply_datasheet(self):
        # get a list of the parameters we will be working with
        sample_count = 0
        meta_cat_count = 0
        for i in range(len(self.dss_objects)):
            changed = False
            dss_obj = self.dss_objects[i]
            print(f'Processing {dss_obj.name}')
            ser = self.sample_meta_info_df.loc[dss_obj.name]
            for meta_cat in self.meta_category_titles:
                print(f'\tChecking {meta_cat}')
                if pd.notnull(ser.at[meta_cat]):
                    # We need to do a float test if we are checking the lat and lon rather than a standard (string)
                    # equality test
                    if meta_cat == 'collection_latitude' or meta_cat == 'collection_longitude':
                        if float(ser.at[meta_cat])!= float(999):
                            if float(getattr(dss_obj, meta_cat)) != float(ser.at[meta_cat]):
                                changed = True
                                print(f'\tChanging {meta_cat} from '
                                      f'{float(getattr(dss_obj, meta_cat))} to {float(ser.at[meta_cat])}')
                                setattr(dss_obj, meta_cat, float(ser.at[meta_cat]))
                                meta_cat_count += 1
                            else:
                                print(f'\t{meta_cat} is already {float(ser.at[meta_cat])}; No change')
                        else:
                            print(f'\tProvided {meta_cat} is invalid for {dss_obj.name} and has been set to {float(999)}; No change')
                    else:
                        if ser.at[meta_cat] != 'NoData':
                            if getattr(dss_obj, meta_cat) != ser.at[meta_cat]:
                                changed = True
                                print(f'\tChanging {meta_cat} from {getattr(dss_obj, meta_cat)} to {ser.at[meta_cat]}')
                                setattr(dss_obj, meta_cat, ser.at[meta_cat])
                                meta_cat_count += 1
                            else:
                                print(f'\t{meta_cat} is already {ser.at[meta_cat]}; No change')
                        else:
                            print(f'\tProvided {meta_cat} is invalid for {dss_obj.name} and has been set to {"NoData"}; No change')
                else:
                    print(f'\t{meta_cat} is null')
            if changed:
                print(f'\tSaving changes to {dss_obj.name}')
                dss_obj.save()
                sample_count += 1
        print('Application complete')
        print(f'Changed {meta_cat_count} meta information category instances across {sample_count} samples')

    def _make_data_sheet_df(self):
        """This will create a curate the df that will hold the information that is to be applied to the
        DataSetSamples of the DataSet in question.
        Much of this code will be similar to the _create_and_check_datasheet method of data_loading.py.
        However, there will be some differences. For example if any of the provided values have to be
        converted back to 999 then we will return an error and not perform the mapping as there
        is no value from applying 999 or null data to samples in the database. We will also allow there to be
        empty columns. The only column that will be mandatory will be the sample_names columns.
        Before applying anychanges to any of the database objects we should ensure that all of these quality
        controls are complete."""
        # Create a pandas df from the data_sheet provided
        # Allow the data_sheet to be in a .csv format or .xlsx format.
        # The sample_meta_df that is created from the data_sheet should be identical irrespective of whether a .csv
        # or a .xlsx is submitted.
        # 1 Check that each of the sample names are unique and that they match an existing sample in
        # the DataSet in question
        # List the sample names that do not match
        # If lat and lon are provided. Check that the lat lon file exists and is in a format that can be used. If not convert to 999.
        # Check that each of the other values can be converted to str

        # Datasheet to pandas df
        self._df_from_path()

        # check sample names are unique
        self._check_vals_of_col_unique(column_name='sample_name')

        self.sample_meta_info_df['sample_name'] = self.sample_meta_info_df['sample_name'].astype(str)
        self.sample_meta_info_df['sample_name'] = self.sample_meta_info_df[
            'sample_name'].str.rstrip().str.lstrip().str.replace(' ', '_').str.replace('/', '_')

        self.sample_meta_info_df.set_index('sample_name', inplace=True, drop=True)

        # ensure sample names are strings
        self.sample_meta_info_df.index = self.sample_meta_info_df.index.map(str)

        self._check_for_binomial()

        self._replace_null_vals_in_meta_info_df()

        self._check_lat_long()

        self._check_vars_can_be_string()

    def _check_vars_can_be_string(self):
        """First convert each of the columns to type string.
        Then make sure that all of the vals are genuine vals of NoData
        """
        for col in ['sample_type', 'host_phylum', 'host_class', 'host_order', 'host_family', 'host_genus',
                        'host_species', 'collection_depth', 'collection_date']:
            self.sample_meta_info_df[col] = self.sample_meta_info_df[col].astype(str)
        for i, sample_name in enumerate(self.sample_meta_info_df.index.values.tolist()):
            for col in ['sample_type', 'host_phylum', 'host_class', 'host_order', 'host_family', 'host_genus',
                        'host_species', 'collection_depth', 'collection_date']:
                try:
                    value = str(self.sample_meta_info_df.at[sample_name, col])
                    if value == 'nan':
                        self.sample_meta_info_df.at[sample_name, col] = 'NoData'
                except:
                    self.sample_meta_info_df.at[sample_name, col] = 'NoData'

    def _set_lat_lon_to_999(self, sample_name):
        self.sample_meta_info_df.at[sample_name, 'collection_latitude'] = float(999)
        self.sample_meta_info_df.at[sample_name, 'collection_longitude'] = float(999)

    def _check_lat_long(self):
        # check the lat long value for each sample listed
        for i, sample_name in enumerate(self.sample_meta_info_df.index):
            dirty_lat = self.sample_meta_info_df.at[sample_name, 'collection_latitude']
            dirty_lon = self.sample_meta_info_df.at[sample_name, 'collection_longitude']
            try:
                clean_lat, clean_lon = check_lat_lon(lat=dirty_lat, lon=dirty_lon)
            except RuntimeError:
                print(
                    f'Unable to convert the Lat Lon values of {sample_name} to decimal degrees.'
                    f'Values will be set to 999'
                )
                print(f'No changes will be made to DataSetSample object {sample_name} for lat or lon\n')
                self._set_lat_lon_to_999(sample_name)
                continue

            # final check to make sure that the values are in a sensible range
            if (-90 <= clean_lat <= 90) and (-180 <= clean_lon <= 180):
                self.sample_meta_info_df.at[sample_name, 'collection_latitude'] = clean_lat
                self.sample_meta_info_df.at[sample_name, 'collection_longitude'] = clean_lon
            else:
                print(f'The lat and lon values for {sample_name} are not in a sensible range '
                      f'({clean_lat}, {clean_lon}). Values will be set to 999')
                print(f'No changes will be made to DataSetSample object {sample_name} for lat or lon\n')
                self._set_lat_lon_to_999(sample_name)
        # finally make sure that the lat and long cols are typed as float
        self.sample_meta_info_df['collection_latitude'] = self.sample_meta_info_df['collection_latitude'].astype(float)
        self.sample_meta_info_df['collection_longitude'] = self.sample_meta_info_df['collection_longitude'].astype(float)

    def _replace_null_vals_in_meta_info_df(self):
        self.sample_meta_info_df = self.sample_meta_info_df.replace('N/A', NaN)\
            .replace('NA', NaN).replace('na', NaN).replace('n/a', NaN)

    def _check_for_binomial(self):
        """People were putting the full binomial in the speices colums. This crops this back to just the
        species component of binomial"""
        for row_name in self.sample_meta_info_df.index.values.tolist():
            current_species_val = self.sample_meta_info_df.at[row_name, 'host_species']
            if not pd.isnull(current_species_val):
                if ' ' in current_species_val:
                    new_species_val = current_species_val.split(' ')[-1]
                    print(f'changing {current_species_val} to {new_species_val} for {row_name}')
                    self.sample_meta_info_df.at[row_name, 'host_species'] = new_species_val

    def _df_from_path(self):
        if self.datasheet_path.endswith('.xlsx'):
            self.sample_meta_info_df = pd.read_excel(
                io=self.datasheet_path, header=0, usecols='A:N', skiprows=[0])
        elif self.datasheet_path.endswith('.csv'):
            with open(self.datasheet_path, 'r') as f:
                data_sheet_as_file = [line.rstrip() for line in f]
            if data_sheet_as_file[0].split(',')[0] == 'sample_name':
                self.sample_meta_info_df = pd.read_csv(
                    filepath_or_buffer=self.datasheet_path)
            else:
                self.sample_meta_info_df = pd.read_csv(
                    filepath_or_buffer=self.datasheet_path, skiprows=[0])
        else:
            sys.exit('Data sheet: {} is in an unrecognised format. '
                     'Please ensure that it is either in .xlsx or .csv format.')

    def _check_vals_of_col_unique(self, column_name):
        # check to see that the values held in a column are unique
        sample_name_counter = Counter(self.sample_meta_info_df[column_name].values.tolist())
        non_unique_name_list = []
        for col_val, count in sample_name_counter.items():
            if count != 1:
                non_unique_name_list.append(col_val)
        if non_unique_name_list:
            print(f'There appear to be non unique {column_name}s in your datasheet:')
            for col_val in non_unique_name_list:
                print(f'\t{col_val}')
            print('Please rectify this in your datasheet and reattempt loading')
            sys.exit()

    def _apply_datasheet_to_dataset_samples(self):
        foo = 'bar'

class CreateStudyAndAssociateUsers:
    """
    Currently there are two situations in which we may want to be creating a new Study object
    As part of this we may also want to be creating new User objects to associate to it
    The first instance is when we are doing a loading operation.
    The second is when we are updating citations
    This Class provides the user input series to allow us to create a new Study object
    and optionally create Users to associate to it.
    """
    def __init__(
            self, ds=None, citation=None, date_time_str=None, list_of_dss_objects=None,
            is_cron_loading=None, study_name=None, study_user_string=None
    ):
        # We should be provided either a DataSet or a citation
        # if ds provided then we are doing a loading operation
        # else if citation provided then we are updating citations
        if ds is not None and citation is not None:
            raise RuntimeError('Please provide either a DataSet object or a citation argument')
        if ds is not None:
            self.dataset_object = ds
            self.citation = None
            assert(date_time_str is not None)
            self.date_time_str = date_time_str
            assert(list_of_dss_objects is not None)
            self.list_of_dss_objects = list_of_dss_objects
            self.is_cron_loading = is_cron_loading
            self.study_name = study_name
            self.study_user_string = study_user_string
        else:
            self.dataset_object = None
            self.citation = citation
            self.date_time_str = str(datetime.utcnow()).split('.')[0].replace('-','').replace(' ','T').replace(':','')
            self.list_of_dss_objects = None
        self.restart = False
        self.study = None

    def create_study_and_user_objects(self):
        """
        When is_cron_loading is True, this is a loading initiated by one of the cron jobs
        We will bypass the user input queries and use the provided name a user string when is_cron_loading
        """
        while True:
            restart = False
            print("\nYou can enter '4' at any of the command prompts to start this process again")
            if self.dataset_object is not None:
                if self.is_cron_loading:
                    study_name = self.study_name
                    continue_text = 'y'
                else:
                    continue_text = input(f'Associate DataSet {self.dataset_object.name} to a Study [y/n]: ')
                    study_name = self.dataset_object.name
            else:
                # If we are working with citation then we already know that the user wants to create a new
                # Study object
                while True:
                    print('The citation is:')
                    print(f'\ttitle: {self.citation["title"]}')
                    print(f'\tauthors: {self.citation["authors"]}')
                    print(f'\tyear: {self.citation["year"]}')
                    study_name = input('Please provide a name for the new Study: ')
                    continue_text = input(f'Is {study_name} correct? [y/n]: ')
                    if continue_text == 'y':
                        break
            if continue_text == 'y':
                if not self._study_of_name_exists(study_name=study_name):
                    new_study_object, restart = self._create_study_object(name=study_name, restart=restart)
                else:
                    new_study_object, restart = self._create_study_object_name_exists(restart, study_name)

                # If we were given a '4', restart the process
                if restart:
                    continue

                while True:
                    if self.is_cron_loading:
                        continue_text = 'y'
                    else:
                        continue_text = input('Associate one or more users to this Study? [y/n]: ')
                    if continue_text == 'y':
                        if self.is_cron_loading:
                            continue_text = self.study_user_string
                        else:
                            continue_text = input("Enter a comma separated string of the usernames for the users you "
                                                  "wish to create. For example: 'jbloggs,vshnazzy': ")
                        user_names_to_associate = continue_text.split(',')
                        if len(user_names_to_associate) == 0:
                            # Then it was left blank
                            # Go back to the prompt
                            print('Please enter at least one username.')
                            continue
                        users_that_already_exist, users_that_do_not_exist = self._check_to_see_if_users_exist(
                            user_names_to_associate)
                        self._print_out_users(users_that_already_exist, users_that_do_not_exist)
                        if self.is_cron_loading:
                            continue_text = 'y'
                        else:
                            continue_text = input("\nIs this correct? [y/n]: ")
                        if continue_text == 'y':
                            self._create_study_users_and_associate(new_study_object, users_that_already_exist,
                                                                   users_that_do_not_exist)
                            break
                        elif continue_text == 'n':
                            # Then something went wrong and we should start again
                            continue
                    elif continue_text == 'n':
                        self.study = self._create_study_associate_dss_objects(new_study_object=new_study_object)
                        print(f'Created {new_study_object} with no users associated.')
                        break
                    elif continue_text == '4':
                        restart = True
                        break

                if restart:
                    continue
                else:
                    # If we are not restarting then we have successfuly associated any necessary
                    # Study and User objects and we can return
                    if self.dataset_object is not None:
                        print(f'\n\nDataSet {self.dataset_object.name} now associated to Study {new_study_object}')
                    else:
                        print(f'citation now associated to Study {new_study_object}')
                    print(f'Study {new_study_object} now associated to user(s):')
                    for user in new_study_object.user_set.all():
                        print(f'\t{user}')
                    return
            elif continue_text == 'n':
                return
            elif continue_text == '4':
                continue
            else:
                print('Unrecognised response. Please answer y or n.')

    def _create_study_users_and_associate(self, new_study_object, users_that_already_exist, users_that_do_not_exist):
        # Now we can create the Study
        self.study = self._create_study_associate_dss_objects(new_study_object)

        self._create_and_associate_users(new_study_object, users_that_already_exist, users_that_do_not_exist)

    def _create_and_associate_users(self, new_study_object, users_that_already_exist, users_that_do_not_exist):
        # First make the new users
        if users_that_do_not_exist:
            print('Creating new users:')
            for user in users_that_do_not_exist:
                new_user = User(name=user)
                new_user.save()
                new_user.studies.add(new_study_object)
                new_user.save()
                print(f'\t{new_user}')
        # Then associate the study to existing Users
        if users_that_already_exist:
            for user in users_that_already_exist:
                user.studies.add(new_study_object)
                user.save()
        # At this point we have completed the associations

    def _create_study_associate_dss_objects(self, new_study_object):
        new_study_object.save()
        if self.dataset_object is not None:
            new_study_object.data_set_samples.set(self.list_of_dss_objects)
            new_study_object.save()
        else:
            while True:
                associate_data_set = input('Associate a DataSet to the study? [y/n]: ')
                if associate_data_set == 'y':
                    uids = input('Please input a comma separated string of the dataset UIDs. E.g. 43,56: ')
                    try:
                        ds_to_associate = DataSet.objects.filter(id__in=[int(_) for _ in uids.split(',')])
                        # Associate the assocaited datasetsamples to the Study
                        self.list_of_dss_objects = DataSetSample.objects.filter(
                            data_submission_from__in=ds_to_associate
                        )
                        new_study_object.data_set_samples.set(self.list_of_dss_objects)
                        new_study_object.save()
                    except Exception as e:
                        print(e)
                        print('\nOne or more of the datasets could not be found')
                        continue
                    new_study_object.save()
                    break
                elif associate_data_set == 'n':
                    break
        return new_study_object
    def _print_out_users(self, users_that_already_exist, users_that_do_not_exist):
        if users_that_already_exist:
            print('The following users already exist and the Study will be associated to them:')
            for user in users_that_already_exist:
                print(f'\t{user}')
        if users_that_do_not_exist:
            print('\nThe following users will be created:')
            for user in users_that_do_not_exist:
                print(f'\t{user}')

    def _check_to_see_if_users_exist(self, list_of_usernames):
        # This will be the actual database ORM User objects
        exists = []
        # This will be a list of usernames to create as strings
        do_not_exist = []
        for username in list_of_usernames:
            try:
                exists.append(User.objects.get(name=username))
            except ObjectDoesNotExist:
                do_not_exist.append(username)
        return exists, do_not_exist

    def _create_study_object_name_exists(self, restart, study_name):
        new_study_object = None
        while True:
            # If we are doing the loading then give the option to overwrite the current
            # DataSetSamples with the new DataSetSamples of the loaded DataSet
            # Otherwise, if we are doing citations, ask if the user wants to associate this citation
            # to the Study. This will mean changing the title etc. to match the citation title
            print(f"Study with name '{study_name}' already exists.")
            if self.dataset_object is not None:
                if self.is_cron_loading:
                    continue_text = 'y'
                else:
                    continue_text = input('Do you want to update/overwrite the current data_set_samples '
                                          'associated with this study to the current DataSetSamples? [y/n]: ')
            else:  # Working with citations
                continue_text = input('Do you want to associate the current citation to this Study?\n'
                                      'This will overwrite the attributes such as title '
                                      'and author_list_string of the '
                                      'Study so that they match those of the citation. [y/n]: ')

            if continue_text == 'y':
                # We will simply set new_study_object to point to the current Study object
                # The new DataSetSample objects will be updated later.
                new_study_object = Study.objects.get(name=study_name)
                break
            elif continue_text == '4':
                restart = True
                break
            elif continue_text == 'n':
                study_name = input('Please provide a name for the new dataset: ')
                if self._study_of_name_exists(study_name):
                    continue
                else:
                    new_study_object, restart = self._create_study_object(restart=restart, name=study_name)
                    break
        return new_study_object, restart

    def _create_study_object(self, restart, name=None):
        while True:
            print('The default Study object that will be created will be:\n')
            if self.dataset_object is not None:
                new_study_object = self._create_study_object_from_dataset_object(name)
            else: # create from the citation object
                new_study_object = self._create_study_object_from_citation(name)
            self._print_non_init_study_summary(new_study_object)
            if self.is_cron_loading:
                continue_text = 'y'
            else:
                continue_text = input('Continue? [y/n]: ')
            if continue_text == 'y':
                break
            elif continue_text == 'n':
                restart = True
                break
            elif continue_text == '4':
                restart = True
                break
            else:
                pass
        return new_study_object, restart

    def _print_non_init_study_summary(self, study):
        print(f'< Study >')
        print(f'name: {study.name}')
        print(f'title: {study.title}')
        print(f'is_published: {study.is_published}')
        print(f'location: {study.location}')
        print(f'run_type: {study.run_type}')
        print(f'article_url: {study.article_url}')
        print(f'data_url: {study.data_url}')
        print(f'data_explorer: {study.data_explorer}')
        print(f'display_online: {study.display_online}')
        print(f'analysis: {study.analysis}')
        print(f'author_list_string: {study.author_list_string}')
        print(f'additional_markers: {study.additional_markers}')
        print(f'creation_time_stamp: {study.creation_time_stamp}')

    def _create_study_object_from_citation(self, name):
        print('Please supply values for the following attributes (leave blank for null):')
        attribute_dict = {}
        attribute_dict['location'] = input('location: ')
        attribute_dict['run_type'] = input('run_type: ')
        attribute_dict['data_url'] = input('data_url: ')
        attribute_dict['additional_markers'] = input('additional_markers: ')
        new_study_object = Study(
            name=name,
            creation_time_stamp=self.date_time_str,
            title=self.citation['title'],
            is_published=True, display_online=True,
            author_list_string=self.citation['authors'], article_url=self.citation['article_url'])
        for k, v in attribute_dict.items():
            if v != '':
                setattr(new_study_object, k, v)
        return new_study_object

    def _create_study_object_from_dataset_object(self, name):
        # We will set analysis to False and then set True when we are
        # doing the output as part of an analysis.
        if name is None:
            new_study_object = Study(
                name=self.dataset_object.name,
                creation_time_stamp=self.date_time_str, title=self.dataset_object.name, analysis=False)
        else:
            new_study_object = Study(
                name=name,
                creation_time_stamp=self.date_time_str, title=name, analysis=False)
        return new_study_object

    def _study_of_name_exists(self, study_name=None):
        try:
            if study_name == None:
                return Study.objects.get(name=self.dataset_object.name)
            else:
                return Study.objects.get(name=study_name)
        except ObjectDoesNotExist:
            return False