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
from hal.mutations.impl.halTreeMutations import runShellCommand

class BedHistogram:

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
        self.genomeLength = None
        self.rate = None
        self.totalEvents = None
        pass
                 
    # read bed file line by line, storing relevant info in class members
    def loadFile(self, bedPath, binSize = 1,
                 events = BedMutations.defaultEvents,
                 hal = None):
        
        self.bins = defaultdict(int)
        self.binSize = binSize

        bm = BedMutations()
        self.totalEvents = 0.0
        
        for line in bm.scan(bedPath, events):
            d = bm.distance()
            if d is not None:
                bin = d / binSize
                self.bins[bin] += 1
                self.totalEvents += 1.0

        # load a "background rate" (number of events / length)
        # if a hal path was given.  we need the hal path for the genome
        # length, since its clumsy to guess the latter from the bed only
        if hal is not None:
            self.__getGenomeLength(hal, bm.genome)
            self.rate = self.totalEvents / self.genomeLength
            print "total events: %d   len: %d    rate: %f    " % (
                self.totalEvents, self.genomeLength, self.rate)


    # compute the genome length using halStats and summing over the
    # sequence lengths
    def __getGenomeLength(self, halPath, genomeName):
        stats = runShellCommand("halStats %s --sequenceStats %s" % (halPath,
                                                                genomeName))
        lines = stats.split("\n")
        self.genomeLength = 0
        for line in lines[1:]:
            tokens = line.split()
            if len(tokens) >= 2:
                self.genomeLength += int(tokens[1].strip(","))
            else:
                assert len(tokens) == 0
        return self.genomeLength

    # compute rate array:  need to review to make better use of numpy!
    def __rateFn(self, x):
        if self.rate is not None:
            total = 0
            y = np.zeros(len(x))
            for xelem in xrange(len(x)):
                s = int(max(0, x[xelem] - (self.binSize / 2.)))
                for i in range(s, s + self.binSize):                    
                    p = self.rate * math.pow(1. - self.rate, i)
                    y[xelem] += int(p * self.totalEvents)
                    total += y[xelem]
            print "rate total %i" % total
            return y

    # draw the histogram of inter-event distances in log scale
    def writeFigure(self, pdfPath):
        pdf = pltBack.PdfPages(pdfPath)
        fig = plt.figure(figsize=(9.0, 4.0), dpi=300, facecolor='w')
        self.__drawData(fig)
        fig.savefig(pdf, format = 'pdf')
        pdf.close()

    # convert bed data into numpy arrays
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
            temp = 500 / self.binSize
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
      
                    
        if self.rate is not None:
            plotlist.append(plt.plot(x, self.__rateFn(x),
                                     color=self.colorList[6],
                                     linestyle='solid', marker='None', 
                                     markeredgecolor=self.colorList[6],
                                     markeredgewidth=0, linewidth=2.5,
                                     markersize=6.0, alpha=0.9)[0])
            legend += ['Background']

        if not self.linearY:
            temp = x[:1000/ self.binSize]
            plotlist.append(plt.plot(temp, np.power(10, model[0] * temp + model[1]),
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
                        default="\"%s\"" % " ".join(BedMutations.defaultEvents),
                        type=str, help="event tags")
    parser.add_argument("--xlimit", default=None, type=int,
                        help="maximum x value to plot")
    parser.add_argument("--hal", default=None, type=str,
                        help="hal file database used for background rate")
    args = parser.parse_args()

    binSize = args.bin
    events =  args.events.split()
    bh = BedHistogram()
    bh.linearY = False
    bh.xlimit = args.xlimit
    bh.loadFile(args.bed, binSize, events, args.hal)
    bh.writeFigure(args.pdf)

if __name__ == "__main__":
    sys.exit(main())
