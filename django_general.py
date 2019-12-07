"""A method that relies on Django feature. E.g. model classes."""
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
from dbApp.models import DataSet, DataAnalysis, DataSetSample
import pandas as pd
import sys
from collections import Counter
from numpy import NaN


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
        self.sample_meta_info_df['collection_latitude'] = self.sample_meta_info_df['collection_latitude'].astype(str)
        self.sample_meta_info_df['collection_longitude'] = self.sample_meta_info_df['collection_longitude'].astype(str)
        for i, sample_name in enumerate(self.sample_meta_info_df.index.values.tolist()):
            lat = self.sample_meta_info_df.at[sample_name, 'collection_latitude']
            lon = self.sample_meta_info_df.at[sample_name, 'collection_longitude']
            lat = lat.rstrip().lstrip().replace(chr(176), '')
            lon = lon.rstrip().lstrip().replace(chr(176), '')

            # 1 - Check to see if we are dealing with nan values
            # any nan values would initially have been numpy.NaN but after convertsion to string above
            # will be "nan".
            # This may cause an issue if the column is in str format already
            if lat == 'nan' or lon == 'nan':
                print(f'Lat and long are currently nan for {sample_name}. Values will be set to 999')
                print(f'No changes will be made to DataSetSample object {sample_name} for lat or lon\n')
                self._set_lat_lon_to_999(sample_name)
                continue
            else:
                # try to see if they are compatable floats
                try:
                    lat_float = float(lat)
                    lon_float = float(lon)
                except Exception:
                    # see if they are decimal degrees only with the hemisphere anotation of degree sign
                    try:
                        if 'N' in lat:
                            lat_float = float(lat.replace('N', '').replace(chr(176), ''))
                            # lat_float should be positive
                            if lat_float < 0:
                                lat_float = lat_float * -1
                        elif 'S' in lat:
                            lat_float = float(lat.replace('S', '').replace(chr(176), ''))
                            # lat_float should be negative
                            if lat_float > 0:
                                lat_float = lat_float * -1
                        else:
                            # There was not an N or S found in the lat so we should raise error
                            raise RuntimeError
                        if 'E' in lon:
                            lon_float = float(lon.replace('E', '').replace(chr(176), ''))
                            # lon_float should be positive
                            if lon_float < 0:
                                lon_float = lon_float * -1
                        elif 'W' in lon:
                            lon_float = float(lon.replace('W', '').replace(chr(176), ''))
                            # lon_float should be negative
                            if lon_float > 0:
                                lon_float = lon_float * -1
                        else:
                            # There was not an N or S found in the lat so we should raise error
                            raise RuntimeError
                    except:
                        # see if they are in proper dms format
                        try:
                            lat_float = self.dms2dec(lat)
                            lon_float = self.dms2dec(lon)
                        # if all this fails, convert to 999
                        except Exception:
                            print(f'Unable to convert the Lat Lon values of {sample_name} to float. Values will be set to 999')
                            print(f'No changes will be made to DataSetSample object {sample_name} for lat or lon\n')
                            self._set_lat_lon_to_999(sample_name)
                            continue
                # final check to make sure that the values are in a sensible range
                if (-90 <= lat_float <= 90) and (-180 <= lon_float <= 180):
                    self.sample_meta_info_df.at[sample_name, 'collection_latitude'] = lat_float
                    self.sample_meta_info_df.at[sample_name, 'collection_longitude'] = lon_float
                else:
                    print(f'The lat and lon values for {sample_name} are not in a sensible range '
                          f'({lat_float}, {lon_float}). Values will be set to 999')
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