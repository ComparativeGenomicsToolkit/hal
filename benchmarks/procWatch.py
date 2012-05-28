#!/usr/bin/env python2.7
""" 
procWatch
28 May 2012, Memorial Day
dent earl, dearl (a) soe ucsc edu

procWatch is a python script to watch a list of processes
and track their max real memory usage, user time and system 
time. 

Because you have better things to do than watching 'top'

"""
import argparse
import os
import sys
import time

def initOptions():
    parser = argparse.ArgumentParser(description='Monitor processes.')
    parser.add_argument('--processes', type=str,
                        help='Comma separated list of process ids to monitor.')
    parser.add_argument('--delay', type=float, default=5.0,
                        help='Processes will be sampled every __ seconds. default=%(default)s.')
    return parser
def checkOptions(parser, args):
    if args.delay < 0:
        parser.error('Delay must be a positive value, not %f' % args.delay)
    if args.processes is None:
        parser.error('Specify --processes')
    args.pids = args.processes.split(',')
    for i, p in enumerate(args.pids, 0):
        try:
            args.pids[i] = int(p)
        except ValueError:
            parser.error('Invalid process id: %s' % p)
        if args.pids[i] < 1:
            parser.error('Invalid process id: %s' % p)
def readMem(f, pid):
    m = 0
    units = ''
    for line in f:
        line = line.strip()
        if not line.startswith('Rss'):
            continue
        data = line.split()
        m += int(data[1])
        if units == '':
            units = data[2]
        else:
            if units != data[2]:
                raise RuntimeError('Units changed in /proc/%d/smaps from %s to %s' % (pid, units, data[2]))
    return m, units
def getCommandLine(pid):
    try:
        f = open(os.path.join('/proc', str(pid), 'cmdline'))
    except IOError:
        return ''
    s = f.read()
    f.close()
    return s
def poll(pids, memory, times):
    success = True
    success = success & pollMemory(pids, memory)
    success = success & pollTimes(pids, times)
    return success
def pollTimes(pids, times):
    read = False
    tick = float(os.sysconf(os.sysconf_names['SC_CLK_TCK']))
    for p in pids:
        try:
            f = open(os.path.join('/proc', str(p), 'stat'))
            read = True
        except IOError:
            continue
        s = f.read()
        f.close
        data = s.split()
        if p not in times:
            times[p] = {}
        times[p]['utime'] = float(data[13]) / tick
        times[p]['stime'] = float(data[14]) / tick
        times[p]['cutime'] = float(data[15]) / tick
        times[p]['cstime'] = float(data[16]) / tick
    return read
def pollCommand(pids, cmds):
    for p in pids:
        cmds[p] = getCommandLine(p)
def pollMemory(pids, memory):
    read = False
    for p in pids:
        try:
            f = open(os.path.join('/proc', str(p), 'smaps'))
            read = True
        except IOError:
            continue
        if p not in memory:
            memory[p] = {}
            memory[p]['max'] = 0
            memory[p]['units'] = ''
        m, units = readMem(f, p)
        if memory[p]['max'] < m:
            memory[p]['max'] = m
            memory[p]['units'] = units
        f.close()
    return read
def getFormattedTime(p, times, mode):
    if mode == 'user':
        d = float(times[p]['utime'])
    elif mode == 'system':
        d = float(times[p]['stime'])
    else:
        raise RuntimeError
    hours = int(d / 60 / 60)
    d -= int(hours * 60 * 60)
    minutes = int(d / 60)
    d -= int(minutes * 60)
    return '%3d:%02d:%05.2f' % (hours, minutes, d)
def report(cmds, memory, times):
    print '%7s %13s %12s %12s %s' % ('PID', 'MaxMem', 'UserTime', 'SysTime', 'Cmd')
    for p in memory:
        print('%7d %10d %s %10s %10s %s' % 
              (p, memory[p]['max'], memory[p]['units'], getFormattedTime(p, times, 'user'), 
               getFormattedTime(p, times, 'system'), cmds[p]))
def monitor(args):
    memory = {}
    times = {}
    cmds = {}
    pollCommand(args.pids, cmds)
    while True:
        if not poll(args.pids, memory, times):
            break
        time.sleep(args.delay)
    report(cmds, memory, times)
        
def main():
    parser = initOptions()
    args = parser.parse_args()
    checkOptions(parser, args)
    monitor(args)

if __name__ == '__main__':
    main()
