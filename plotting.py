from matplotlib.patches import Rectangle  # Rectangle is used despite it being greyed out in pycharm
from matplotlib.collections import PatchCollection
# https://stackoverflow.com/questions/21784641/installation-issue-with-matplotlib-python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from collections import defaultdict
import pandas as pd
from exceptions import DistanceTypeNotIdentifiedError
import os
import numpy as np
import sys
from datetime import datetime
from general import ThreadSafeGeneral
import json
plt.ioff()

class DistScatterPlotter:
    def __init__(self, csv_path, date_time_str):
        self.output_directory = os.path.dirname(csv_path)
        self.clade = self.output_directory.split('/')[-1]
        self.plotting_df = pd.read_csv(csv_path, sep=',', lineterminator='\n', header=0, index_col=0)
        # Check to see that there are more than two samples with sequences from this clade
        # and explicitly check to see that the PC2 exists
        # if not, then raise a run time exception
        if not len(self.plotting_df.index.values.tolist()) > 2 or 'PC2' not in list(self.plotting_df):
            print(f'Insufficient samples or distance between samples for clade {self.clade}. '
                  f'Distances will not be plotted for this clade.')
            raise RuntimeError()
        self.f, self.ax = plt.subplots(1, 1, figsize=(9, 9))
        self.x_values = self.plotting_df['PC1'].values.tolist()[:-1]
        self.y_values = self.plotting_df['PC2'].values.tolist()[:-1]
        self.fig_output_base = None
        self.date_time_str = date_time_str
        self.output_path_list = []
        if 'braycurtis' in csv_path:
            self.dist_type = 'braycurtis'
        elif 'unifrac' in csv_path:
            self.dist_type = 'unifrac'
        else:
            raise DistanceTypeNotIdentifiedError()
        if '_no_sqrt' in csv_path:
            self.sqrt = False
        elif 'sqrt' in csv_path:
            self.sqrt = True
        else:
            raise DistanceTypeNotIdentifiedError

    def create_base_scatter_plot(self):
        self.ax.scatter(self.x_values, self.y_values, c='black', marker='o')

    def _add_proportion_explained_labels(self):
        self.ax.set_xlabel('PC1; explained = {}'.format('%.3f' % self.plotting_df['PC1'][-1]))
        self.ax.set_ylabel('PC2; explained = {}'.format('%.3f' % self.plotting_df['PC2'][-1]))

    def _add_title(self, title_prefix):
        if self.sqrt:
            self.ax.set_title(f'{title_prefix} {self.clade} {self.dist_type} sqrt')
        else:
            self.ax.set_title(f'{title_prefix} {self.clade} {self.dist_type} no_sqrt')

    def _output_dist_scatter(self):
        plt.tight_layout()
        sys.stdout.write('\rsaving as .svg')
        svg_path = '{}.svg'.format(self.fig_output_base)
        plt.savefig(svg_path)
        png_path = '{}.png'.format(self.fig_output_base)
        sys.stdout.write('\rsaving as .png')
        plt.savefig(png_path, dpi=600)
        sys.stdout.write('\r\nDistance plots output to:')
        sys.stdout.write('\n{}'.format(svg_path))
        sys.stdout.write('\n{}\n'.format(png_path))
        self.output_path_list.extend([svg_path, png_path])


class DistScatterPlotterSamples(DistScatterPlotter):
    def __init__(self, csv_path, date_time_str, labels=True):
        super().__init__(csv_path=csv_path, date_time_str=date_time_str)
        self.labels = labels

    def make_sample_dist_scatter_plot(self):
        self.create_base_scatter_plot()
        self._annotate_plot_with_sample_names()
        self._add_title(title_prefix='between sample distances clade')

        if self.sqrt:
            self.fig_output_base = os.path.join(
                self.output_directory,
                f'{self.date_time_str}_between_sample_distances_clade_{self.clade}_{self.dist_type}_sqrt')
        else:
            self.fig_output_base = os.path.join(
                self.output_directory,
                f'{self.date_time_str}_between_sample_distances_clade_{self.clade}_{self.dist_type}_no_sqrt')
        self._output_dist_scatter()

    def _annotate_plot_with_sample_names(self):
        if self.labels:
            sample_names = self._get_sample_names()
            for i, txt in enumerate(sample_names):
                self.ax.annotate(txt, (self.x_values[i], self.y_values[i]))

    def _get_sample_names(self):
        return self.plotting_df.index.values.tolist()[:-1]


class DistScatterPlotterTypes(DistScatterPlotter):
    def __init__(self, csv_path, date_time_str):
        super().__init__(csv_path=csv_path, date_time_str=date_time_str)

    def make_type_dist_scatter_plot(self):
        self.create_base_scatter_plot()
        self._annotate_plot_with_type_uids()
        self._add_title(title_prefix='between its2 type profile distances clade')
        if self.sqrt:
            self.fig_output_base = os.path.join(
                self.output_directory,
                f'{self.date_time_str}_between_its2_type_prof_dist_clade_{self.clade}_{self.dist_type}_sqrt')
        else:
            self.fig_output_base = os.path.join(
                self.output_directory,
                f'{self.date_time_str}_between_its2_type_prof_dist_clade_{self.clade}_{self.dist_type}_no_sqrt')
        self._output_dist_scatter()

    def _annotate_plot_with_type_uids(self):
        for i, txt in enumerate(self.plotting_df.index.values.tolist()[:-1]):
            self.ax.annotate(txt, (self.x_values[i], self.y_values[i]))


class SubPlotter:
    """A class that can be used for the sub plots of the SeqStackedBarPlotter and the TypeStackedBarPlotter.
    """
    def __init__(self, parent_plotter_instance, index_of_this_subplot):
        self.parent_plotter = parent_plotter_instance
        self.patches_list = []
        self.x_tick_label_list = []
        self.x_index_for_plot = 0
        self.colour_list = []
        self.index_of_this_subplot = index_of_this_subplot
        self.subplot_axes = self.parent_plotter.axarr[self.index_of_this_subplot]
        self.end_slice = self._get_end_index_for_slicing_plotting_data()
        self.num_samples_in_this_subplot = len(
                self.parent_plotter.output_count_table_as_df.index.values.tolist()[
                self.index_of_this_subplot * self.parent_plotter.samples_per_subplot:self.end_slice])
        # Make a custom colour map
        # https://matplotlib.org/api/_as_gen/matplotlib.colors.ListedColormap.html
        self.listed_colour_map = None
        self.patches_collection = None

    def plot_seq_subplot(self):
        self._create_rect_patches_and_populate_colour_list()

        self._make_listed_colour_map()

        self._make_patches_collection()

        self._draw_patches_on_axes()

        self._format_axes()

    def _format_axes(self):
        # make it so that the x axes is constant length that will be the num of samples per subplot
        self.subplot_axes.set_xlim(0 - 0.5, self.parent_plotter.samples_per_subplot - 0.5)
        self.subplot_axes.set_ylim(0, 1)
        self.subplot_axes.set_xticks(range(self.num_samples_in_this_subplot))
        self.subplot_axes.set_xticklabels(self.x_tick_label_list, rotation='vertical', fontsize=6)
        self.subplot_axes.spines['right'].set_visible(False)
        self.subplot_axes.spines['top'].set_visible(False)
        # as well as getting rid of the top and right axis splines
        # I'd also like to restrict the bottom spine to where there are samples plotted but also
        # maintain the width of the samples
        # I think the easiest way to do this is to hack a bit by setting the x axis spines to invisible
        # and then drawing on a line at y = 0 between the smallest and largest ind (+- 0.5)
        self.subplot_axes.spines['bottom'].set_visible(False)
        self.subplot_axes.add_line(
            Line2D((0 - 0.5, self.num_samples_in_this_subplot - 0.5), (0, 0), linewidth=2, color='black'))

    def _draw_patches_on_axes(self):
        self.subplot_axes.add_collection(self.patches_collection)
        self.subplot_axes.autoscale_view()
        self.subplot_axes.figure.canvas.draw()

    def _make_patches_collection(self):
        self.patches_collection = PatchCollection(self.patches_list, cmap=self.listed_colour_map)
        self.patches_collection.set_array(np.arange(len(self.patches_list)))

    def _make_listed_colour_map(self):
        self.listed_colour_map = ListedColormap(self.colour_list)

    def _get_end_index_for_slicing_plotting_data(self):
        if self.index_of_this_subplot == self.parent_plotter.number_of_subplots - 1:
            end_slice = self.parent_plotter.num_samples
        else:
            end_slice = self.parent_plotter.samples_per_subplot * (self.index_of_this_subplot + 1)
        return end_slice

    def _create_rect_patches_and_populate_colour_list(self):
        for sample_uid in self.parent_plotter.output_count_table_as_df.index.values.tolist()[
                      self.index_of_this_subplot * self.parent_plotter.samples_per_subplot:self.end_slice]:

            sys.stdout.write(f'\rPlotting sample: {self.parent_plotter.smp_uid_to_smp_name_dict[int(sample_uid)]}')

            self._add_sample_names_to_tick_label_list(sample_uid)
            # for each sample we will start at 0 for the y and then add the height of each bar to this
            bottom = 0
            # for each sequence, create a rect patch
            # the rect will be 1 in width and centered about the ind value.
            non_zero_indices = self.parent_plotter.output_count_table_as_df.loc[sample_uid].to_numpy().nonzero()[0]
            current_sample_series = self.parent_plotter.output_count_table_as_df.loc[sample_uid]
            non_zero_sample_series = current_sample_series.iloc[non_zero_indices]
            sample_total = non_zero_sample_series.sum()
            for ser_index, rel_abund in non_zero_sample_series.iteritems():

                self.patches_list.append(Rectangle(
                    (self.x_index_for_plot - 0.5, bottom),
                    1,
                    rel_abund / sample_total, color=self.parent_plotter.colour_dict[ser_index]))
                self.colour_list.append(self.parent_plotter.colour_dict[ser_index])
                bottom += rel_abund / sample_total

            self.x_index_for_plot += 1

    def _add_sample_names_to_tick_label_list(self, sample_uid):
        sample_name = self.parent_plotter.smp_uid_to_smp_name_dict[int(sample_uid)]
        if len(sample_name) < 20:
            self.x_tick_label_list.append(sample_name)
        else:
            self.x_tick_label_list.append(f'uid_{int(sample_uid)}')


class LegendPlotter:
    """This class can be used by the SeqStackedBarPlotter and the TypeStackedBarPlotter to handle
    the plotting of the legend subplot.
    """
    def __init__(self, parent_plotter, type_plotting=False):
        # whether we are plotting types
        self.type_plotting=type_plotting
        self.parent_plotter = parent_plotter
        self.ax_to_plot_on = self.parent_plotter.axarr[-1]
        # legend setup parameters
        self.y_coord_increments = 100 / self.parent_plotter.max_n_rows
        self.leg_box_depth = 2 / 3 * self.y_coord_increments

        self.x_coord_increments = 100 / self.parent_plotter.max_n_cols
        self.leg_box_width = self.x_coord_increments / 3
        self._set_n_rows_and_last_row_len()
        self.column_count = 0

    def plot_legend_seqs(self):
        self._set_ylim_and_x_lim_and_invert_y_axis()

        self._plot_legend_rows()

        self._remove_frames_from_axis()

    def _plot_legend_rows(self):
        if not self.type_plotting:
            sys.stdout.write(
                f'\nGenerating figure legend for {str(self.parent_plotter.num_leg_cells)} most common sequences\n')
        else:
            sys.stdout.write(
                f'\nGenerating figure legend for {str(self.parent_plotter.num_leg_cells)} most common ITS2 '
                f'type profiles\n')

        for row_increment in range(min(self.n_rows, self.parent_plotter.max_n_rows)):

            if self._this_is_last_row_of_legend(row_increment=row_increment):
                for col_increment in range(self.parent_plotter.max_n_cols):
                    self._plot_legend_row(row_increment=row_increment, col_increment=col_increment)

                    self.column_count += 1
            else:
                for col_increment in range(self.last_row_len):
                    self._plot_legend_row(row_increment=row_increment, col_increment=col_increment)

                    self.column_count += 1

    def _set_ylim_and_x_lim_and_invert_y_axis(self):
        # Once we know the number of rows, we can also adjust the y axis limits
        self.ax_to_plot_on.set_xlim(0, 100)
        self.ax_to_plot_on.set_ylim(0, ((self.n_rows - 1) * self.y_coord_increments) + self.leg_box_depth)
        self.ax_to_plot_on.invert_yaxis()

    def _set_n_rows_and_last_row_len(self):
        if not self.type_plotting:  # we are plotting sequences
            col_elements_to_plot = len(self.parent_plotter.ordered_list_of_seqs_names)
        else:  # we are plotting types
            col_elements_to_plot = len(self.parent_plotter.sorted_type_prof_uids_by_local_abund)
        if col_elements_to_plot < self.parent_plotter.num_leg_cells:
            if col_elements_to_plot % self.parent_plotter.max_n_cols != 0:
                self.n_rows = int(col_elements_to_plot / self.parent_plotter.max_n_cols) + 1
                self.last_row_len = col_elements_to_plot % self.parent_plotter.max_n_cols
            else:
                self.n_rows = int(col_elements_to_plot / self.parent_plotter.max_n_cols)
                self.last_row_len = self.parent_plotter.max_n_cols
        else:
            self.n_rows = self.parent_plotter.max_n_rows
            self.last_row_len = self.parent_plotter.max_n_cols

    def _this_is_last_row_of_legend(self, row_increment):
        return (row_increment + 1) != self.n_rows

    def _remove_frames_from_axis(self):
        self.ax_to_plot_on.set_frame_on(False)
        self.ax_to_plot_on.get_xaxis().set_visible(False)
        self.ax_to_plot_on.get_yaxis().set_visible(False)

    def _plot_legend_row(self, row_increment, col_increment):
        leg_box_x, leg_box_y = self._add_legend_rect(col_increment=col_increment, row_increment=row_increment)
        self._add_legend_text(leg_box_x, leg_box_y)

    def _add_legend_text(self, leg_box_x, leg_box_y):
        text_x = leg_box_x + self.leg_box_width + (0.2 * self.leg_box_width)
        text_y = leg_box_y + (0.5 * self.leg_box_depth)
        if not self.type_plotting:
            self.ax_to_plot_on.text(
                text_x, text_y, self.parent_plotter.ordered_list_of_seqs_names[self.column_count],
                verticalalignment='center', fontsize=8)
        else:
            type_name_to_print = self.parent_plotter.type_uid_to_type_name_dict[
                self.parent_plotter.sorted_type_prof_uids_by_local_abund[self.column_count]]
            if len(type_name_to_print) > 15:
                type_name_to_print = f'{type_name_to_print[:14]}...'

            self.ax_to_plot_on.text(
                text_x, text_y, type_name_to_print,
                verticalalignment='center', fontsize=8)

    def _add_legend_rect(self, col_increment, row_increment):
        leg_box_x = col_increment * self.x_coord_increments
        leg_box_y = row_increment * self.y_coord_increments
        if not self.type_plotting:
            self.ax_to_plot_on.add_patch(Rectangle(
                (leg_box_x, leg_box_y), width=self.leg_box_width, height=self.leg_box_depth,
                color=self.parent_plotter.colour_dict[
                    self.parent_plotter.ordered_list_of_seqs_names[self.column_count]]))
        else:
            self.ax_to_plot_on.add_patch(Rectangle(
                (leg_box_x, leg_box_y), width=self.leg_box_width, height=self.leg_box_depth,
                color=self.parent_plotter.colour_dict[
                    self.parent_plotter.sorted_type_prof_uids_by_local_abund[self.column_count]]))
        return leg_box_x, leg_box_y


class TypeStackedBarPlotter:
    """Class for plotting the type count table output"""
    def __init__(self, type_relative_abund_count_table_path, output_directory, date_time_str=None):
        self.type_rel_abund_count_table_path = type_relative_abund_count_table_path
        self.output_directory = output_directory
        if date_time_str:
            self.date_time_str = date_time_str
        else:
            self.date_time_str = str(datetime.utcnow()).split('.')[0].replace('-','').replace(' ','T').replace(':','')
        self.fig_output_base = os.path.join(self.output_directory, 'its2_type_profiles', f'{self.date_time_str}')
        self.max_n_cols = 5
        self.max_n_rows = 10
        self.num_leg_cells = self.max_n_cols * self.max_n_rows

        self.sorted_type_prof_uids_by_local_abund = None
        self.smp_uid_to_smp_name_dict = None
        self.type_uid_to_type_name_dict = None
        self.output_count_table_as_df = self._create_output_df_and_populate_smpl_id_to_smp_name_dict()
        self.colour_dict = self._set_colour_dict()
        self.num_samples = len(self.output_count_table_as_df.index.values.tolist())
        self.samples_per_subplot = 50
        self.number_of_subplots = self._infer_num_subplots()
        self.f, self.axarr = plt.subplots(self.number_of_subplots + 1, 1, figsize=(10, 3 * self.number_of_subplots))
        self.output_path_list = []

    def plot_stacked_bar_profiles(self):
        print('\n\nPlotting ITS2 type profile abundances')
        for sub_plot_index in range(self.number_of_subplots):
            sub_plotter = SubPlotter(index_of_this_subplot=sub_plot_index, parent_plotter_instance=self)
            sub_plotter.plot_seq_subplot()

        self._plot_legend()

        plt.tight_layout()

        self._write_out_plot()


    def _write_out_plot(self):
        sys.stdout.write('\nFigure generation complete')
        sys.stdout.write('\nFigures output to:\n')

        svg_path = f'{self.fig_output_base}_type_abundance_stacked_bar_plot.svg'
        self.output_path_list.append(svg_path)
        sys.stdout.write(f'{svg_path}\n')
        plt.savefig(svg_path)

        png_path = f'{self.fig_output_base}_type_abundance_stacked_bar_plot.png'
        self.output_path_list.append(png_path)
        sys.stdout.write(f'{png_path}\n')
        plt.savefig(png_path, dpi=600)

    def _plot_legend(self):
        legend_plotter = LegendPlotter(parent_plotter=self, type_plotting=True)
        legend_plotter.plot_legend_seqs()

    def _infer_num_subplots(self):
        # number of subplots will be one per smp_per_plot
        # and if tehre are remainers be sure to add an extra plot for this
        if (self.num_samples % self.samples_per_subplot) != 0:
            return int(self.num_samples / self.samples_per_subplot) + 1
        else:
            return int(self.num_samples / self.samples_per_subplot)

    def _create_output_df_and_populate_smpl_id_to_smp_name_dict(self):
        sp_output_df = pd.read_csv(
            self.type_rel_abund_count_table_path, sep='\t', lineterminator='\n', skiprows=[1, 2, 3, 5], header=None)

        # get a list of tups that are the profile uid and the abundances zipped together
        type_profile_to_abund_tup_list = [
            (name, int(abund)) for name, abund in
            zip(sp_output_df.iloc[0][2:].values.tolist(), sp_output_df.iloc[1][2:].values.tolist())]

        # convert the names that are numbers into int strings rather than float strings.
        int_temp_list = []
        for name_abund_tup in type_profile_to_abund_tup_list:
            try:
                int_temp_list.append((str(int(name_abund_tup[0])), int(name_abund_tup[1])))
            except:
                int_temp_list.append((name_abund_tup[0], int(name_abund_tup[1])))
        type_profile_to_abund_tup_list = int_temp_list

        # need to drop the rows that contain the sequence accession and species descriptions
        index_to_drop_from = None
        for i, row_name in enumerate(sp_output_df.iloc[:, 0]):
            if 'Sequence accession' in str(row_name):
                # then we want to drop all rows from here until the end
                index_to_drop_from = i
                break

        sp_output_df = sp_output_df.iloc[:index_to_drop_from]

        # now make a dict of sample id to sample name so that we can work with uids
        self.smp_uid_to_smp_name_dict = {int(smp_uid): str(smp_name) for smp_uid, smp_name in zip(sp_output_df.iloc[3:, 0], sp_output_df.iloc[3:, 1])}

        # now make a dict of of type id to type name so that we can also work eith uids for the types
        # as there could be types with identical names
        self.type_uid_to_type_name_dict = {int(type_uid): type_name for type_uid, type_name in zip(sp_output_df.iloc[0, 2:], sp_output_df.iloc[2, 2:])}

        # now drop the sample name columns
        sp_output_df.drop(columns=1, inplace=True)

        # now drop the type name columns
        sp_output_df.drop(index=2, inplace=True)

        # now promote the its2_type_prof names to columns headers and drop the local abund row and .
        sp_output_df.columns = ['sample_id'] + [int(a[0]) for a in type_profile_to_abund_tup_list]
        sp_output_df.drop(index=[0, 1], inplace=True)

        self.sorted_type_prof_uids_by_local_abund = [
            int(a[0]) for a in sorted(type_profile_to_abund_tup_list, key=lambda x: x[1], reverse=True)]

        # convert the sample_id col to numeric
        sp_output_df['sample_id'] = pd.to_numeric(sp_output_df['sample_id'])

        return sp_output_df.set_index(keys='sample_id', drop=True).astype('float')


    def _set_colour_dict(self):
        """The colour dict was created when making the profile rectagle array for the
        javascript output. We need to maintain syncronisity so we will read this in.
        JSON automatically sets key values to strings so we will need to convert the strings back to ints"""
        with open(os.path.join(self.output_directory, 'html', 'prof_color_dict.json'), 'r') as f:
            str_color_dict = json.load(fp=f)

        return {int(k):v for k, v in str_color_dict.items()}


class SeqStackedBarPlotter():
    """Class for plotting the sequence count table output"""
    def __init__(
            self, seq_relative_abund_count_table_path_post_med, seq_relative_abund_df_pre_med, output_directory,
            no_pre_med_seqs, ordered_seq_list, date_time_str=None, ordered_sample_uid_list=None):
        self.seq_relative_abund_count_table_path_post_med = seq_relative_abund_count_table_path_post_med
        self.seq_relative_abund_df_pre_med = seq_relative_abund_df_pre_med
        self.ordered_seq_list = ordered_seq_list
        self.root_output_directory = output_directory
        self.post_med_output_directory = os.path.join(self.root_output_directory, 'post_med_seqs')
        self.pre_med_output_directory = os.path.join(self.root_output_directory, 'pre_med_seqs')
        if date_time_str:
            self.date_time_str = date_time_str
        else:
            self.date_time_str = str(datetime.utcnow()).split('.')[0].replace('-','').replace(' ','T').replace(':','')
        self.fig_output_base = os.path.join(self.post_med_output_directory, f'{self.date_time_str}')
        self.smp_uid_to_smp_name_dict = None
        self.output_count_table_as_df = pd.read_csv(
            self.seq_relative_abund_count_table_path_post_med,
            sep='\t', lineterminator='\n', header=0, index_col=0
        )
        self._format_output_df_and_populate_smpl_id_to_smp_name_dict()
        self.ordered_list_of_seqs_names = self._set_ordered_list_of_seqs_names()
        # legend parameters and vars
        self.max_n_cols = 8
        self.max_n_rows = 7
        self.num_leg_cells = self.max_n_rows * self.max_n_cols
        self.thread_safe_general = ThreadSafeGeneral()
        self.colour_dict = self.thread_safe_general.set_seq_colour_dict(self.ordered_list_of_seqs_names)
        # plotting vars
        self.ordered_sample_uid_list = self._set_ordered_sample_uid_list_and_reorder_df(ordered_sample_uid_list)
        self.num_samples = len(self.output_count_table_as_df.index.values.tolist())
        self.samples_per_subplot = 50
        self.number_of_subplots = self._infer_number_of_subplots()
        # we add  1 to the n_subplots here for the legend at the bottom
        self.f, self.axarr = plt.subplots(self.number_of_subplots + 1, 1, figsize=(10, 3 * self.number_of_subplots))
        self.output_path_list = []
        self.no_pre_med_seqs = no_pre_med_seqs

    def plot_stacked_bar_seqs(self):
        self._plot_stacked_bar_seqs_post_med()
        if not self.no_pre_med_seqs:
            if self.seq_relative_abund_df_pre_med is not None: # This will be None if the premed seq output was not done.
                self._plot_stacked_bar_seqs_pre_med()

    def _plot_stacked_bar_seqs_post_med(self):
        print('\n\nPlotting sequence abundances')
        for sub_plot_index in range(self.number_of_subplots):
            sub_plotter = SubPlotter(index_of_this_subplot=sub_plot_index, parent_plotter_instance=self)
            sub_plotter.plot_seq_subplot()

        self._plot_legend()

        plt.tight_layout()

        self._write_out_plot()

        self.output_path_list.extend(
            [
                f'{self.fig_output_base}_seq_abundance_stacked_bar_plots.svg',
                f'{self.fig_output_base}_seq_abundance_stacked_bar_plots.png'
            ])

    def _write_out_plot(self):
        sys.stdout.write('\nFigure generation complete')
        sys.stdout.write('\nFigures output to:\n')

        svg_path = f'{self.fig_output_base}_seq_abundance_stacked_bar_plot.svg'
        sys.stdout.write(f'{svg_path}\n')
        plt.savefig(svg_path)

        png_path = f'{self.fig_output_base}_seq_abundance_stacked_bar_plot.png'
        sys.stdout.write(f'{png_path}\n')
        plt.savefig(png_path, dpi=600)

    def _plot_legend(self):
        legend_plotter = LegendPlotter(parent_plotter=self)
        legend_plotter.plot_legend_seqs()

    def _get_end_index_for_slicing_plotting_data(self, i):
        if i == self.number_of_subplots - 1:
            end_slice = self.num_samples
        else:
            end_slice = self.samples_per_subplot * (i + 1)
        return end_slice

    def _add_sample_names_to_tick_label_list(self, sample, x_tick_label_list):
        sample_name = self.smp_uid_to_smp_name_dict[int(sample)]
        if len(sample_name) < 20:
            x_tick_label_list.append(self.smp_uid_to_smp_name_dict[int(sample)])
        else:
            x_tick_label_list.append(f'uid_{int(sample)}')

    def _infer_number_of_subplots(self):
        if (self.num_samples % self.samples_per_subplot) != 0:
            number_of_subplots = int(self.num_samples / self.samples_per_subplot) + 1
        else:
            number_of_subplots = int(self.num_samples / self.samples_per_subplot)
        return number_of_subplots

    def _set_ordered_sample_uid_list_and_reorder_df(self, ordered_sample_uid_list):
        """If we are plotting this in companion with an ITS2 type profile output then we will be passed a
        ordered_sample_uid_list. It is very useful to have the ITS2 type profile output figure and the seq figure
        in the same sample order for direct comparison.
        If this output is not associated with an ITS2 type profile output then we will need to
        generate the sample order from scratch.
        In theory the output should already be somewhat ordered in that the samples should be in order of similarity.
        However, these have the artifical clade ordering (sorted by clade and then by abundance of seqs) so for the
        plotting it will probably be better to get a new
        order for the samples that is not constrained to the order of the clades. For this we should order as usual
        according to the most common majority sequences and then within this grouping we should order according to the
        the abundance of these sequences within the samples.
        """
        if not ordered_sample_uid_list:
            self.ordered_sample_uid_list = self._generate_sample_order_de_novo()
        else:
            self.ordered_sample_uid_list = ordered_sample_uid_list

        self._reorder_df_by_new_sample_and_seq_order()

        return self.ordered_sample_uid_list

    def _reorder_df_by_new_sample_and_seq_order(self):
        self.output_count_table_as_df = self.output_count_table_as_df[self.ordered_list_of_seqs_names]
        self.output_count_table_as_df = self.output_count_table_as_df.reindex(
            [int(a) for a in self.ordered_sample_uid_list])

    def _generate_sample_order_de_novo(self):
        """At this stage we have the ordered list of seqs we now need to order the samples
        this method will return us the names of the samples in order that they should be plotted.
        """
        # {sequence_name_found_to_be_most_abund_in_sample: num_samples_it_was_found_to_be_most_abund_in}
        max_seq_ddict = defaultdict(int)
        # {most_abundant_seq_name: [(dss.id, rel_abund_of_most_abund_seq) for samples with that seq as most abund]}
        seq_to_samp_ddict = defaultdict(list)

        # for each sample get the columns name of the max value of a div not including the columns in the following:
        no_maj_seq = []
        for sample_id_to_sort in self.output_count_table_as_df.index.values.tolist():
            smp_series = self.output_count_table_as_df.loc[sample_id_to_sort].astype('float')

            rel_abund_of_max_abund_seq_name = smp_series.max()
            if not rel_abund_of_max_abund_seq_name > 0:
                no_maj_seq.append(sample_id_to_sort)
            else:
                max_abund_seq_name = smp_series.idxmax()
                # add a tup of sample name and rel abund of seq to the seq_to_samp_dict
                seq_to_samp_ddict[max_abund_seq_name].append((sample_id_to_sort, rel_abund_of_max_abund_seq_name))
                # add this to the ddict count
                max_seq_ddict[max_abund_seq_name] += 1

        # then once we have compelted this for all sequences
        # generate the sample order according to the sequence order
        ordered_sample_list = []
        # get an ordered list of the sequencs according to the max_seq_ddict
        ordered_list_of_sequences = [x[0] for x in sorted(max_seq_ddict.items(), key=lambda x: x[1], reverse=True)]

        for seq_to_order_samples_by in ordered_list_of_sequences:
            tup_list_of_samples_that_had_sequence_as_most_abund = seq_to_samp_ddict[seq_to_order_samples_by]
            ordered_list_of_samples_for_seq_ordered = [
                x[0] for x in sorted(tup_list_of_samples_that_had_sequence_as_most_abund,
                                     key=lambda x: x[1], reverse=True)]
            ordered_sample_list.extend(ordered_list_of_samples_for_seq_ordered)

        ordered_sample_list.extend(no_maj_seq)

        return ordered_sample_list

    def _set_ordered_list_of_seqs_names(self):
        """Get a list of the sequences in order of their abundance and use this list to create the colour dict
        The abundances are got by simply summing up the columns
        """
        abundance_dict = {}
        for col in list(self.output_count_table_as_df):
            abundance_dict[col] = sum(self.output_count_table_as_df[col])

        # get the names of the sequences sorted according to their totalled abundance
        return [x[0] for x in sorted(abundance_dict.items(), key=lambda x: x[1], reverse=True)]

    def _format_output_df_and_populate_smpl_id_to_smp_name_dict(self):
        """Drop the QC columns from the SP output df and also drop the clade summation columns
        we will be left with just columns for each one of the sequences found in the samples
        we need to drop the rows first before we can make the smp_id_to_smp_name_dict else
        we will have the final row names in the index which are not convertable to int
        need to make the smp_id_to_smp_name_dict before dropping the sample_name col"""

        self._drop_meta_info_rows()

        # sample names to string
        self.output_count_table_as_df['sample_name'] = self.output_count_table_as_df['sample_name'].astype(str)

        self._populate_smpl_id_to_smp_name_dict()

        self._drop_non_seq_abund_cols_and_set_df_types()

    def _drop_non_seq_abund_cols_and_set_df_types(self):
        self.output_count_table_as_df.drop(
            columns=[col for col in list(self.output_count_table_as_df) if col not in self.ordered_seq_list]
            , inplace=True
        )
        self.output_count_table_as_df = self.output_count_table_as_df.astype('float')
        self.output_count_table_as_df.index = self.output_count_table_as_df.index.astype('int')

    def _populate_smpl_id_to_smp_name_dict(self):
        self.smp_uid_to_smp_name_dict = {
            int(uid): smp_name for uid, smp_name in
            zip(self.output_count_table_as_df.index.values.tolist(), self.output_count_table_as_df['sample_name'].values.tolist())}

    def _drop_meta_info_rows_from_df(self, meta_index_to_cut_from):
        self.output_count_table_as_df.drop(index=self.output_count_table_as_df.index[range(meta_index_to_cut_from, 0, 1)], inplace=True)

    def _drop_meta_info_rows(self):
        # In order to be able to drop the DIV row at the end and the meta information rows, we should
        # drop all rows that are after the DIV column. We will pass in an index value to the .drop
        # that is called here. To do this we need to work out which index we are working with
        meta_index_to_cut_from = None
        index_values_as_list = self.output_count_table_as_df.index.values.tolist()
        for i in range(-1, -(len(index_values_as_list)), -1):
            if index_values_as_list[i].startswith('seq'):
                # then this is the index (in negative notation) that we need to cut from
                meta_index_to_cut_from = i
                break
        self.output_count_table_as_df.drop(index=self.output_count_table_as_df.index[range(meta_index_to_cut_from, 0, 1)], inplace=True)

    @staticmethod
    def _get_colour_list():
        colour_list = [
            "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#FFDBE5", "#7A4900",
            "#0000A6",
            "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400",
            "#4FC601",
            "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
            "#B903AA",
            "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
            "#372101",
            "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
            "#0CBD66",
            "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
            "#BEC459",
            "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
            "#FF913F",
            "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7",
            "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
            "#201625",
            "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
            "#CB7E98",
            "#A4E804", "#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489",
            "#806C66",
            "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5",
            "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58", "#7A7BFF", "#D68E01",
            "#353339",
            "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A",
            "#001325",
            "#02525F", "#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75",
            "#8D8546",
            "#9695C5", "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
            "#00005F",
            "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01", "#6B94AA", "#51A058",
            "#A45B02",
            "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406",
            "#F4D749",
            "#0045D2", "#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE",
            "#C6DC99",
            "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
            "#C6005A",
            "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183", "#AA5199", "#B5D6C3",
            "#A38469",
            "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03",
            "#E4FFFC",
            "#1BE177", "#BCB1E5", "#76912F", "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213",
            "#A76F42",
            "#89412E", "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
            "#BDC9D2",
            "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00", "#061203", "#DFFB71",
            "#868E7E",
            "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C",
            "#00B57F",
            "#545C46", "#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"]
        return colour_list

    def _plot_stacked_bar_seqs_pre_med(self):
        pre_med_seq_plotter = self.PreMedSeqPlotter(parent=self)
        pre_med_seq_plotter.plot_stacked_bar_seqs()

    class PreMedSeqPlotter():
        def __init__(self, parent):
            self.parent = parent
            self.thread_safe_general = ThreadSafeGeneral()
            self.root_output_directory = self.parent.root_output_directory
            self.fig_output_base = os.path.join(self.parent.pre_med_output_directory, f'{self.parent.date_time_str}')
            self.smp_uid_to_smp_name_dict = None
            self.output_count_table_as_df = self._curate_output_count_table(self.parent.seq_relative_abund_df_pre_med)
            self.ordered_list_of_seqs_names = list(self.output_count_table_as_df)
            # legend parameters and vars
            self.max_n_cols = 8
            self.max_n_rows = 7
            self.num_leg_cells = self.max_n_rows * self.max_n_cols
            self.colour_dict = self.thread_safe_general.set_seq_colour_dict_w_reference_c_dict(self.ordered_list_of_seqs_names, self.parent.colour_dict)
            # plotting vars
            self.num_samples = len(self.output_count_table_as_df.index.values.tolist())
            self.samples_per_subplot = 50
            self.number_of_subplots = self._infer_number_of_subplots()
            # we add  1 to the n_subplots here for the legend at the bottom
            self.f, self.axarr = plt.subplots(self.number_of_subplots + 1, 1, figsize=(10, 3 * self.number_of_subplots))
            self.output_path_list = []

        def _curate_output_count_table(self, rel_abund_df):
            self.smp_uid_to_smp_name_dict = {
                uid: str(name) for uid, name in
                zip(rel_abund_df.index.values.tolist(), rel_abund_df['sample_name'].values.tolist())}
            df = rel_abund_df.drop(columns='sample_name')
            df = df.astype('float')
            df = df.reindex(
                [int(a) for a in self.parent.ordered_sample_uid_list])
            return df.astype('float')

        def plot_stacked_bar_seqs(self):
            print('\n\nPlotting sequence abundances')
            for sub_plot_index in range(self.number_of_subplots):
                sub_plotter = SubPlotter(index_of_this_subplot=sub_plot_index, parent_plotter_instance=self)
                sub_plotter.plot_seq_subplot()

            self._plot_legend()

            plt.tight_layout()

            self._write_out_plot()

            self.output_path_list.extend(
                [
                    f'{self.fig_output_base}_seq_abundance_stacked_bar_plots.svg',
                    f'{self.fig_output_base}_seq_abundance_stacked_bar_plots.png'
                ])

        def _write_out_plot(self):
            sys.stdout.write('\nFigure generation complete')
            sys.stdout.write('\nFigures output to:\n')

            svg_path = f'{self.fig_output_base}_seq_abundance_stacked_bar_plot.svg'
            sys.stdout.write(f'{svg_path}\n')
            plt.savefig(svg_path)

            png_path = f'{self.fig_output_base}_seq_abundance_stacked_bar_plot.png'
            sys.stdout.write(f'{png_path}\n')
            plt.savefig(png_path, dpi=600)

        def _plot_legend(self):
            legend_plotter = LegendPlotter(parent_plotter=self)
            legend_plotter.plot_legend_seqs()

        def _infer_number_of_subplots(self):
            if (self.num_samples % self.samples_per_subplot) != 0:
                number_of_subplots = int(self.num_samples / self.samples_per_subplot) + 1
            else:
                number_of_subplots = int(self.num_samples / self.samples_per_subplot)
            return number_of_subplots