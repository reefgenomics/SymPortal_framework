import os
import shutil
from multiprocessing import Queue, Manager, Process
import sys

from django import db
class SPDataAnalysis:
    def __init__(self, workflow_manager_parent, data_analysis_obj):
        self.parent = workflow_manager_parent
        self.temp_wkd = os.path.join(self.parent.symportal_root_directory)
        self._setup_temp_wkd()
        self.data_analysis_obj = data_analysis_obj
        # The abundance that a given DIV must be found at when it has been considered 'unlocked'
        # https://github.com/didillysquat/SymPortal_framework/wiki/The-SymPortal-logic#type-profile-assignment---logic
        self.unlocked_abundance = 0.0001
        self.clade_list = list('ABCDEFGH')
        self.ccs_of_analysis = self.data_analysis_obj.get_clade_collections()
        # List that will hold a dictionary for each clade
        # Each dictionary will hold key = footprint (set of sequences)
        # value = [[] []] where []0 = list of cladeCollections containing given footprint
        # and []1 = list of majority sequence for given sample
        self.clade_footp_dicts_list = [{} for _ in self.clade_list]



    def analyse_data(self):
        print('Beginning profile discovery')
        footprint_dict_pop_handler = FootprintDictPopHandler(sp_data_analysis_parent=self)
        footprint_dict_pop_handler.populate_clade_footprint_dicts()


    def _setup_temp_wkd(self):
        if os.path.exists(self.temp_wkd):
            shutil.rmtree(self.temp_wkd)
        os.makedirs(self.temp_wkd, exist_ok=True)

class FootprintDictPopHandler:
    """Will handle the execusion of the FootprintDictWorker. This worker will populate the
    SPDataAnalysis master_cladal_list_of_footpinrt_dicts"""
    def __init__(self, sp_data_analysis_parent):
        self.parent = sp_data_analysis_parent
        self.cc_mp_queue = Queue()
        self.output_mp_queue = Queue()
        self._populate_queue()
        self.all_procs = []

    def _populate_queue(self):
        for cc in self.parent.ccs_of_analysis:
            self.cc_mp_queue.put(cc)

        for N in range(self.parent.parent.args.num_proc):
            self.cc_mp_queue.put('STOP')

    def populate_clade_footprint_dicts(self):

        # close all connections to the db so that they are automatically recreated for each process
        # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
        db.connections.close_all()

        for n in range(self.parent.parent.args.num_proc):
            p = Process(target=self._start_footprint_dict_workers,
                        args=())
            self.all_procs.append(p)
            p.start()

        self._collect_cc_info()

    def _collect_cc_info(self):
        kill_number = 0
        while 1:
            cc_info_holder = self.output_mp_queue.get()
            if cc_info_holder == 'kill':
                kill_number += 1
                if kill_number == self.parent.parent.args.num_procs:
                    break
            else:
                # TODO we got here and then the bubkin turned up!
                if cc_info_holder.footprint not in self.parent.clade_footp_dicts_list[
                    cc_info_holder.clade_index]:
                    self.parent.clade_footp_dicts_list[
                        cc_info_holder.clade_index][
                        cc_info_holder.footprint] = FootprintRepresentative(
                        cc=cc_info_holder.cc, cc_maj_ref_seq=cc_info_holder.maj_ref_seq)
                else:
                    self.parent.clade_footp_dicts_list[cc_info_holder.clade_index][cc_info_holder.footprint].cc_list.append(
                        cc_info_holder.cc)
                    self.parent.clade_footp_dicts_list[cc_info_holder.clade_index][
                        cc_info_holder.footprint].maj_seq_list.append(
                        cc_info_holder.maj_ref_seq)

        # First wait for the workers to finish
        for p in self.all_procs:
            p.join()

    def _start_footprint_dict_workers(self):
        for cc in iter(self.cc_mp_queue.get, 'STOP'):

            footprint = cc.cutoff_footprint(self.parent.parent.data_analysis_obj.within_clade_cutoff)

            self.output_mp_queue.put(CCInfoHolder(
                footprint=footprint,
                clade_index=self.parent.parent.clade_list.index(cc.clade), cc=cc,
                maj_ref_seq=cc.maj()))

            sys.stdout.write(f'\rFound footprint {footprint}')

        self.output_mp_queue.put('kill')



class FootprintRepresentative:
    def __init__(self, cc, cc_maj_ref_seq):
        self.cc_list = [cc]
        self.maj_seq_list = [cc_maj_ref_seq]

class CCInfoHolder:
    """An object used in the FootprintDictPopWorker to hold information for each of the clade collections
    that will be used to popualte the clade footprint dicts"""
    def __init__(self, footprint, clade_index, cc, maj_ref_seq):
        self.footprint = footprint
        self.clade_index = clade_index
        self.cc = cc
        self.maj_ref_seq = maj_ref_seq