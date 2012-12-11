#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Make histogram of inter-event distances in BED file
"""
import argparse
import os
import sys
import copy
import random
import math
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator, LogFormatter # minor tick marks
import matplotlib.mlab as mlab


from hal.analysis.neutralIndel.bedMutations import BedMutations

class BedHistogram:

    def __init__(self):
        self.linearY = True
        pass
                 
    # read bed file line by line, storing relevant info in class members
    def loadFile(self, bedPath, binSize = 1,
                 events = BedMutations.defaultEvents):
        
        self.bins = defaultdict(int)
        self.binSize = binSize

        bm = BedMutations()
        
        for line in bm.scan(bedPath, events):
            d = bm.distance()
            if d is not None:
                bin = d / binSize
                self.bins[bin] += 1

    def writeFigure(self, pdfPath):
        pdf = pltBack.PdfPages(pdfPath)
        fig = plt.figure(figsize=(9.0, 4.0), dpi=300, facecolor='w')
        self.__drawData(fig)
        fig.savefig(pdf, format = 'pdf')
        pdf.close()
        
    def __extractPlotTables(self):
        keylist = self.bins.keys()
        keylist.sort()
        x = []
        y = []
        for key in keylist:
            x.append(key)
            y.append(self.bins[key])
        # center the x values in the middle of their bin
        x = (np.array(x) * self.binSize) - self.binSize / 2.0
        y = np.array(y)
        return x, y 

    def __drawData(self, fig):
        ax = fig.add_axes([0.1, 0.2, 0.85, 0.7])
        x, y = self.__extractPlotTables()
        
        colorList = ['#1f77b4', # dark blue
                     '#aec7e8', # light blue
                     '#ff7f0e', # bright orange
                     '#ffbb78', # light orange
                     '#4B4C5E', # dark slate gray
                     '#9edae5', # light blue 
                     '#7F80AB', # purple-ish slate blue
                     '#c7c7c7', # light gray
                     '#9467bd', # dark purple
                     '#c5b0d5', # light purple
                     '#d62728', # dark red
                     '#ff9896', # light red
                     ]
        plotlist = []
       
        plotlist.append(plt.plot(x, y, color=colorList[0],
                                 linestyle='none', marker='.', 
                                 markeredgecolor=colorList[0],
                                 markeredgewidth=0, linewidth=0.5,
                                 markersize=10.0, alpha=0.8)[0])
        
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        if xmin < -1.0:
            plt.xlim(0 - .02 * xmax, xmax * 1.02)
        else:
            plt.xlim(xmax - ((xmax - xmin) * 1.02), xmax * 1.02)
        pltTitle = plt.title("Title")
        plt.setp(pltTitle, fontsize='x-small') # legend fontsize
        for loc, spine in ax.spines.iteritems():
            if loc in ['left','bottom']:
                # outward by 10 points
                spine.set_position(('outward', 10)) 
            elif loc in ['right','top']:
                # don't draw spine               
                spine.set_color('none') 
            else:
                raise ValueError('unknown spine location: %s' % loc)
        # turn off ticks where there is no spine
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        if not self.linearY:
            ax.set_yscale('log')
            ax.yaxis.set_minor_locator(LogLocator(base=10, subs=range(1, 10)))
        plt.xlabel('Distance')
        plt.ylabel('Count')

def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    parser.add_argument("bed", help="input bed")
    parser.add_argument("pdf", help="output pdf")
    parser.add_argument("--bin", default=10, type=int, help="bin size")
    args = parser.parse_args()

    binSize = args.bin
    events =  BedMutations.defaultEvents
    bh = BedHistogram()
    bh.loadFile(args.bed, binSize, events)
    bh.writeFigure(args.pdf)

if __name__ == "__main__":
    sys.exit(main())
