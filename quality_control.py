
from util import argparse
from auxiliary.QualityControl import QualityControl

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', dest='qc_input', nargs='+', default=["/home/tsupinie/hootpy/data/Surface_METAR_20110311_0000.nc"])
    ap.add_argument('--output', dest='qc_output', nargs='+', default=[])
    ap.add_argument('--rules', dest='qc_rules', nargs='+', default=[])

    args = ap.parse_args()

    if args.qc_output == []:
        args.qc_output = None

    qc = QualityControl(args.qc_rules)
    qc.qualityControl(args.qc_input, args.qc_output)

    return

if __name__ == "__main__":
    main()
