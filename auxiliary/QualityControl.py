import numpy as np

from util.util import fatalError, warning
from util.dataIO import DataIO

import re

class QualityControl:
    def __init__(self, rules):
        """
        __init__()
        Purpose:    Constructor for the Quality Control class.  Handles the parsing of the rule strings (parses out the arguments list so it can 
                        determine what variables it needs from the file), and calls eval() to convert to a Python function object.
        Parameters: rules [type=list,tuple]
                        A list of strings that can be evaluated to Python function objects (easiest way is the lambda keyword).
        """
        self.rules = {}

        if len(rules) == 0:
            warning("No rules given for quality control.")
            return

        for rule in rules:
            # Eliminate multiple spaces between "lambda" and the opening parenthesis of the argument list (screws up the next line if there are extra spaces)
            rule = re.sub(r"(?<=lambda)[\s]+", " ", rule)

            # Parse out the argument list to the function
            match = re.search(r"(?<=lambda )[\w, ]+(?= *:)", rule)

            # Get the variables from the match object
            variables = match.group(0)

            # Evaluate the rule (should be a python lambda function) and stick it in the dictionary, with the variables in the argument list as the key
            self.rules[variables] = eval(rule)

        return

    def qualityControl(self, file_names, out_file_names=None):
        """
        qualityControl() [public]
        Purpose:    Main interface for the quality control utility.  Will quality-control a list of files.  If out_file_names is specified, each element 
                        in the list will be matched pairwise with those in file_names.  So the output from the quality control of file_names[0] will be put
                        in out_file_names[0], and likewise for each index.  If it is not specified, each file in file_names will be overwritten with quality
                        controlled data.
        Parameters: file_names [type=list,tuple]
                        A list or tuple of strings giving file names to read in and quality control.
                    out_file_names [type=list,tuple]
                        An optional argument for specifying where to write output from the quality control to.
        Returns:    [nothing]
        """
        if type(file_names) not in [ list, tuple ]:
            file_names = [ file_names ]

        if out_file_names is None:
            out_file_names = [ '' for f in file_names ]

        if len(file_names) != len(out_file_names):
            fatalError("The number of file names must be the same as the number of output file names.")

        for file_name, out_file_name in zip(file_names, out_file_names):
            self._qc(file_name, out_file_name)
        return

    def _qc(self, file_name, out_file_name=''):
        """
        _qc() [private]
        Purpose:    Carries out a quality control operation on a single file.
        Parameters: file_name [type=str]
                        The name of the file to quality control.
                    out_file_name [type=str]
                        Optional file name specifying where to write the output of the quality control procedure.  If not specified, the file given by 
                        file_name will be overwritten with quality-controlled data.
        Returns:    [nothing]
        """

        file = DataIO(file_name=file_name, mode="rw")
        file_data_good = {} 
        file_data_dims = ()

        for variables, function in self.rules.iteritems():
            # Loop through all the rules that we've defined and find the indexes of all the good data.
            if variables.find(",") > -1:
                variables = re.split(", *", variables)
            else:
                variables = [ variables ]

            for var in variables:
                if var not in file_data_good:
                    file_data_dims = file.get_variable_attribute(var, 'dimensions')
                    file_data_good[var] = np.array(np.ones(tuple([ file.get_dimension(d) for d in file_data_dims ]), dtype=bool))

            # Get the variable data from the file (returns as a list of numpy arrays)
            file_data = file.get_variable(variables)

            # Plug the list straight into the quality control function.  Do the element-wise "and" procedure between the file_data_good array and the returned array immediately.
            good = function(*file_data)

            for var in variables:
                file_data_good[var] &= good

#       file_data_good &= file_data_good[:,5].data.reshape((file.get_dimension(file_data_dims[0]), 1))

        for var in file_data_good.keys():

            # Do the splicing for every variable with the same dimensions as the variables we've qc'ed.
            qc_data = np.where(file_data_good[var], file.get_variable(var), np.nan)
            file.set_variable(var, qc_data)

#           new_dims = []
#           keep_indexes = []
#           for idx, var_dim in enumerate(var_dims):
#               if var_dim in file_data_dims:
#                   new_dims.append("%s_qc" % var_dim)
#                   keep_indexes.append(np.unique(np.where(file_data_good)[idx]))
#               else:
#                   new_dims.append(var_dim)
#                   keep_indexes.append(np.arange(file.get_dimension(var_dim), dtype=np.int32))

#           file.splice_variable(variable, tuple(new_dims), tuple(keep_indexes))

        # Sew 'er up
        file.close(out_file_name=out_file_name)
        return

if __name__ == "__main__":
    qc = QualityControl("/home/tsupinie/hootpy/data/Surface_METAR_20110311_0000.nc")
