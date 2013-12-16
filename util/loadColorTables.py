"""
loadColorTables.py
Purpose:  Utility script that loads GEMPAK color table files and converts them to HootPy fconfig files.
Started:  David John Gagne on March 5, 2012
"""

def main():
    import argparse
    parser = argparse.ArgumentParser(prog='loadColorTables',description='Load a GEMPAK or matplotlib color table and convert it into an fconfig file.')
    parser.add_argument('cname',metavar='COLORNAME',
                        help='Matplotlib colormap name or GEMPAK colortable file name')
    parser.add_argument('--type',choices=['M','G'],default='M',
                        help='Specify whether colortable is from Matplotlib (M) or GEMPAK (G).')
    parser.add_argument('-s','--start',type=float,required=True,
                        help='Starting value of interval')
    parser.add_argument('-e','--end',type=float,required=True,
                        help='Ending value of interval')
    parser.add_argument('-i','--interval',type=float,default=0,
                        help='Standard interval between contours.  Do not include if you wish plot a specified number of contours')
    parser.add_argument('-n','--ncontours',type=int,default=0,
                        help='Number of contours.  Do not include if you want to have a standard interval.')
    parser.add_argument('-l','--sample',action='store_true',
                        help='Sample colors evenly across the map.')
    parser.add_argument('-o','--offset',type=int,default=0,
                        help='Colormap offset value. Use if you want your color range to start at a different value')
    parser.add_argument('-x','--extend',choices=['both','max','min','neither'],default='both',
                        help='Specify whether you want to extend the colormap past the max, min, or both extrema')
    parser.add_argument('-w','--wrap',action='store_true',
                        help='If specified, the colortable will wrap around to the beginning+offset')
    parser.add_argument('-u','--unit',default='',
                        help='Units of contours')
    parser.add_argument('-f','--file',default=None,
                        help='Output file name')
    parser.add_argument('-r','--reverse',action='store_true',
                        help='Reverse colormap')
    args = parser.parse_args()
    fill = makeFConfigFill(args.cname,args.type,args.start,args.end,interval=args.interval,num_contours=args.ncontours,offset=args.offset,extend=args.extend,wrap=args.wrap,reverse=args.reverse,sample=args.sample)
    if args.file is not None:
        writeFConfigFile(args.file,args.unit,fill)
    else:
        for f in fill:
            print f

def loadGempakColorTable(filename):
    """loadGempakColorTable
       Purpose:  Read a GEMPAK colortable file and convert the colors to a list of hex strings
       Parameters:  filename
       Returns:  list of hex colors
    """
    import struct
    ct_file = open(filename)
    color_tuples = []
    hex_colors = []
    for line in ct_file:
        if '!' in line[0]:
            continue
        else:
            color_line = line.split()
            color_start = [x.isdigit() for x in color_line].index(True)
            color_tuples.append(tuple([int(x) for x in color_line[color_start:color_start+3]]))
            hex_colors.append('#' + struct.pack('BBB',*color_tuples[-1]).encode('hex'))
    return hex_colors

def loadMPColorTable(cname):
    import struct
    from matplotlib.colors import rgb2hex
    from pylab import get_cmap
    cmap = get_cmap(cname)
    if cmap is None:
        print 'Error: colormap %s does not exist' % cname
        exit()
    c_range = range(0,256)
    hex_colors = []
    for c in c_range:
        hex_colors.append(rgb2hex(cmap(c)))
    return hex_colors

def makeFConfigFill(cname,type,start,end,interval=0,num_contours=0,offset=0,extend='both',wrap=True,reverse=False,sample=False):
    """makeFConfigFill()
       Purpose:  Load the color table file and convert it into the HootPy fill list
       Parameters:
            cname: Either GEMPAK colortable file name or matplotlib colortable name
            type:  'G' for GEMPAK color table or 'M' for matplotlib color table
            start: first fill value
            end:  last fill value
            interval:  if > 0, then the distance between values will be specified by the interval.  Otherwise
                       if < 1, then the values between start and end will be evenly spaced over the length of
                           the colortable.
            extend:  (Default 'both'):  Can be 'both','min','max','neither'.  Specifies whether the beginning and end colors
                     should be used for all values outside the specified range or if they should be ignored.
            wrap:  (Default True):  If True, then the color table is reused from the beginning if more values are required.
       Returns:  fill - a list of fill interval dictionaries describing the color and the upper and lower limits of its range.
    """
    import numpy as np
    fill = []
    if type.upper()=='G':
        hex_colors = loadGempakColorTable(cname)
    elif type.upper()=='M':
        hex_colors = loadMPColorTable(cname)
    else:
        print 'Error:  Improper type specified'
        exit()
    print hex_colors
    if reverse:
        hex_colors.reverse()
    if interval > 0:
        values = np.arange(start,end+interval,interval)
    elif num_contours > 0:
        values = np.linspace(start,end,num_contours)
    else:
        values = np.linspace(start,end,len(hex_colors))
    c = offset
    if sample:
        c_interval = len(hex_colors) / len(values)
    else:
        c_interval = 1
    for i,v in enumerate(values):
        if (extend=='both' or extend=='min') and i==0:
            fill.append({'lower_bound':None,'upper_bound':np.round(values[i],2),'color':hex_colors[c]})
            c+= c_interval
        elif (extend=='both' or extend=='max') and i==len(values) - 1:
            fill.append({'lower_bound':np.round(values[i],2),'upper_bound':None,'color':hex_colors[c]})
        if i < len(values) - 1:
            fill.append({'lower_bound':np.round(values[i],2),'upper_bound':np.round(values[i+1],2),'color':hex_colors[c]})
        if c >= 0 and c < len(hex_colors):
            c+= c_interval
        if wrap and c == len(hex_colors):
            c = 0
        elif not wrap and c == len(hex_colors):
            c = -1
    
    return fill

def writeFConfigFile(config_filename,unit,fill):
    fconfig_file = open(config_filename,'w')
    fconfig_file.write('#Auto-generated fconfig file\n')
    unit_str = 'units = "%s"\n' % unit
    fconfig_file.write(unit_str)
    fconfig_file.write('\n')
    fconfig_file.write('fill = [\n')
    for row in fill[:-1]:
        fconfig_file.write('\t' + str(row) + ',\n')
    fconfig_file.write('\t' + str(fill[-1]) + '\n')
    fconfig_file.write(']\n')
    fconfig_file.close()

def testMain():
    ct_filename = '/usr/local/nawips/gempak/tables/luts/ir_drgb.tbl'
    hex_colors = loadGempakColorTable(ct_filename)
    fill = makeFConfigFill(ct_filename,50,400,10)
    fill2 = makeFConfigFill(ct_filename,50,400)
    writeFConfigFile('../config/ir_sat.fconfig','W m-2',fill2)
    for f in fill2:
        print f
    print hex_colors
    ct_filename2 = '/usr/local/nawips/gempak/tables/colors/coltbl.xwp'
    hex_colors2 = loadGempakColorTable(ct_filename2)
    print hex_colors2

if __name__=="__main__":
    main()
