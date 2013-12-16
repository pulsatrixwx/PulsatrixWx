#!/usr/bin/python

import util.argparse as argparse
from util.config import Config
#from util.callCount import dumpCallCounts
from base.MapPlotter import MapPlotter

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', dest='config_filename', required=True)
    ap.add_argument('--meta', dest='meta_filename')

    args = ap.parse_args()

    config = Config(args.config_filename)

    if args.meta_filename:
        config['meta_filename'] = args.meta_filename

    map = MapPlotter(config)
    map.loadData()
    map.plot()

#   dumpCallCounts()

    return

if __name__ == "__main__":
    main()
