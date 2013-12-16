
from util import argparse
from util.dataIO import DataIO
import types

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-file', dest='data_file', default="data/ECMWF_Global_2p5_20120203_0000.nc")
    ap.add_argument('--attributes', dest='attrs', action='store_true')
    ap.add_argument('--dimensions', dest='dims', action='store_true')
    ap.add_argument('--variables', dest='vars', action='store_true')
    ap.add_argument('--variable', dest='var', default="")

    args = ap.parse_args()

    data_file = DataIO(file_name=args.data_file)
    if args.dims:
        print "Dimensions:"
        dims = data_file.get_dimensions().items()
        max_length = max([len(n[0]) for n in dims])
        for dim_name, dim in dims:
            print "  %s%s" % (dim_name, " " * (max_length - len(dim_name))), len(dim)

    if args.attrs:
        print "File Attributes:"
        attr_list = data_file.get_file_attribute_list()
        max_length = max([len(n) for n in attr_list])
        for attribute in sorted(attr_list):
            print "  %s%s" % (attribute, " " * (max_length - len(attribute))), data_file.get_file_attribute(attribute)

    if args.vars:
        print "Variables:"
        for variable in sorted(data_file.get_variable_list()):
            print "  %s" % variable

    if args.var != "":
        print "Attributes for Variable %s" % args.var
        attr_list = data_file.get_variable_attribute_list(args.var)
        max_length = max([len(n) for n in attr_list])
        for attribute in sorted(attr_list):
            print "  %s%s" % (attribute, " " * (max_length - len(attribute))), data_file.get_variable_attribute(args.var, attribute)
        dim_list = data_file.get_variable_dimensions(args.var)
        print "Dimensions:  ",dim_list
    return

if __name__ == "__main__":
    main()
