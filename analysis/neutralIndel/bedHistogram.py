#!/usr/bin/env python3

#Copyright (C) 2012 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

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
from hal.mutations.impl.halTreeMutations import runShellCommand
from hal.analysis.neutralIndel.backgroundRate import getBackgroundRate

class BedHistogram(object):

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
    def __init__(self):
        self.linearY = True
        self.xlimit = None
        self.totalEvents = None
        self.bgRate = None
        self.title = "Title"
                 
    # read bed file line by line, storing relevant info in class members
    def loadFile(self, bedPath, binSize = 1, bgRate = None,
                 events = BedMutations.defaultEvents):

        self.bins = defaultdict(int)
        self.binSize = binSize
        self.bgRate = bgRate

        bm = BedMutations()
        self.totalEvents = 0.0
        
        for line in bm.scan(bedPath, events):
            d = bm.distance()
            if d is not None:
                bin = d / binSize
                self.bins[bin] += 1
                self.totalEvents += 1.0

    # compute rate array:  need to review to make better use of numpy!
    def __rateFn(self, x):
        if self.bgRate is not None:
            total = 0
            y = np.zeros(len(x))
            for xelem in range(len(x)):
                s = int(max(0, x[xelem] - (self.binSize / 2.)))
                for i in range(s, s + self.binSize):                    
                    p = self.bgRate * math.pow(1. - self.bgRate, i)
                    y[xelem] += int(p * self.totalEvents)
                    total += y[xelem]
            print("rate total %i" % total)
            return y

    # draw the histogram of inter-event distances in log scale
    def writeFigure(self, pdfPath, title):
        self.title = title
        pdf = pltBack.PdfPages(pdfPath)
        fig = plt.figure(figsize=(9.0, 4.0), dpi=300, facecolor='w')
        self.__drawData(fig)
        fig.savefig(pdf, format = 'pdf')
        pdf.close()

    # convert bed data into numpy arrays
    def __extractPlotTables(self):
        keylist = list(self.bins.keys())
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

    # this is mostly code from Dent 
    def __drawData(self, fig):
        ax = fig.add_axes([0.1, 0.2, 0.85, 0.7])
        x, y = self.__extractPlotTables()

        # truncate along x-axis if desired
        if self.xlimit is not None:
            x = x[:self.xlimit / self.binSize]
            y = y[:self.xlimit / self.binSize]

        # compute linear fit and r2
        if not self.linearY:
            ylog = np.log10(y)
            xcof = np.vstack((x, np.ones(len(x)))).T
            temp = 900 / self.binSize
            model, resid = np.linalg.lstsq(xcof[:temp], ylog[:temp])[:2]
            r2 = 1 - resid / (ylog.size * ylog.var())

        plotlist = []
        legend = []
       
        plotlist.append(plt.plot(x, y, color=self.colorList[1],
                                 linestyle='none', marker='.', 
                                 markeredgecolor=self.colorList[1],
                                 markeredgewidth=0, linewidth=0.5,
                                 markersize=6.0, alpha=0.75)[0])
        legend += ['Data']
      
                    
        if self.bgRate is not None:
            plotlist.append(plt.plot(x, self.__rateFn(x),
                                     color=self.colorList[6],
                                     linestyle='solid', marker='None', 
                                     markeredgecolor=self.colorList[6],
                                     markeredgewidth=0, linewidth=2.5,
                                     markersize=6.0, alpha=0.9)[0])
            legend += ['Background']

        if not self.linearY:
            temp = x[:100000/ self.binSize]
            plotlist.append(plt.plot(temp,
                                     np.power(10, model[0] * temp + model[1]),
                                     color=self.colorList[2],
                                     linestyle='solid', marker='None', 
                                     markeredgecolor=self.colorList[2],
                                     markeredgewidth=0, linewidth=2.5,
                                     markersize=10.0, alpha=0.9)[0])
            legend += ['Linear Fit r2=%s' % r2]
        
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        if xmin < -1.0:
            plt.xlim(0 - .02 * xmax, xmax * 1.02)
        else:
            plt.xlim(xmax - ((xmax - xmin) * 1.02), xmax * 1.02)
        pltTitle = plt.title(self.title)
        plt.setp(pltTitle, fontsize='x-small') # legend fontsize
        for loc, spine in ax.spines.items():
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
            ax.yaxis.set_minor_locator(LogLocator(base=10, subs=list(range(1, 10))))
            ax.yaxis.set_data_interval(0.01, None)
        plt.xlabel('Distance')
        plt.ylabel('Count')
        plt.ylim(ymin = 0.5)

        leg = plt.legend(plotlist, legend, loc=1, numpoints=1)
        plt.setp(leg.get_texts(), fontsize='x-small') # legend fontsize
        leg._drawFrame = False
        
    
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    parser.add_argument("bed", help="input bed")
    parser.add_argument("pdf", help="output pdf")
    parser.add_argument("--bin", default=10, type=int, help="bin size")
    parser.add_argument("--events",
                        default=" ".join(BedMutations.defaultEvents),
                        type=str, help="event tags")
    parser.add_argument("--xlimit", default=None, type=int,
                        help="maximum x value to plot")
    parser.add_argument("--backgroundBed", default=None, type=float,
                        help="regions for computing background rate")
    parser.add_argument("--title", default=None, type=str,
                        help="title of plot")
    args = parser.parse_args()

    binSize = args.bin
    events =  args.events.split()
    bh = BedHistogram()
    bh.linearY = False
    bh.xlimit = args.xlimit

    bgRate = None
    if args.backgroundBed is not None:
        count, size = getBackgroundRate(args.bed, args.backgroundBed, events)
        bgRate = float(count) / float(size)
    
    bh.loadFile(args.bed, binSize, bgRate, events)
    bh.writeFigure(args.pdf, args.title)

if __name__ == "__main__":
    sys.exit(main())
