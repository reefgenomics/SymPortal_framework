#!/usr/bin/env python3

# from dbApp.models import ReferenceSequence
import sys
import json

class SeqMatcher:
    def __init__(self):
        # Dict: k= nucleotide sequence, v=dict (k=DataSetSample object, key=absolute abund of seq in dss)
        with open(sys.argv[1], 'r') as f:
            self.seq_list = json.load(f)
        # Dict: k= nucleotide sequence, v=corresponding ReferenceSequence object
        with open(sys.argv[2], 'r') as f:
            self.rs_list = json.load(f)
        # Make a list for faster parseing
        self.rs_set = set(self.rs_list)
        # The full path to which the match and non-match dicts should be output via compress pickle
        self.match_dict_output_path = sys.argv[3]
        self.non_match_list_output_path = self.match_dict_output_path.replace('match_dict', 'non_match_list')
        # The dictionary to be returned k=nucleotide sequence, v=ReferenceSequence object
        self.match_dict = {}
        self.non_match_list = []

    def match(self):
        for nuc_seq in self.seq_list:
            if not self._match_found(nuc_seq):
                self.non_match_list.append(nuc_seq)

        with open(self.match_dict_output_path, 'w') as f:
            json.dump(self.match_dict, f)
        with open(self.non_match_list_output_path, 'w') as f:
            json.dump(self.non_match_list, f)

    def _match_found(self, nuc_seq):
        # Try to match the exact sequence
        if nuc_seq in self.rs_set:
            self.match_dict[nuc_seq] = nuc_seq
            return True
        else:
            # Finally try to find a super or sub match
            for rs_seq in self.rs_list:
                if nuc_seq in rs_seq or rs_seq in nuc_seq:
                    # Then this is a match
                    self.match_dict[nuc_seq] = rs_seq
                    return True
        return False

if __name__ == '__main__':
    SeqMatcher().match()
