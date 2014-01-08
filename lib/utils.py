#!/usr/bin/env python


import sys
import os
import fnmatch
import re
import os
import pwd
import collections


def get_FilesList(path,pattern=None):
    """
    traverse a dir recursively and find files with full paths
    that match the given pattern
    """
    
    if pattern is None:
        print '[%s]: No pattern supplied' % (SUB)
        return
    
    SUB='get_FilesList'
    found_files = []
    for root,dirs,files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                filename = str(os.path.join(root,basename))
                found_files.append(filename.strip())
    print '[%s]: Found %d files at %s' % ( SUB,len(found_files),path)
    return found_files





def _get_logger(logfile,dual_mode=True,logger_name=None):
    """
        return a logger object which
        which writes to both consol and logfile
        
        Reference:
        got from http://docs.python.org/2/howto/logging-cookbook.html (Logging to multiple destinations)
    """
        
    import logging
    
    if logger_name is None:
        logger_name == __name__

    
    logging.basicConfig(filename=logfile,
                        level=logging.DEBUG,
                        filemode='w',
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m/%d/%Y %H:%M')
    
    
    
    #define a Handler which writes INFO messages or higher to the sys.stderr
    console_output = logging.StreamHandler()
    console_output.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console_output.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console_output)
    
    logger = logging.getLogger(logger_name)
    
    return logger





def chunks(list,n):
    '''
        yield successive n-sized chunks from list
    '''
    for i in xrange(0,len(list),n):
        yield list[i:i+n]


def uniqify_list(list):
   # Not order preserving
   keys = {}
   for e in list:
       keys[e] = 1
   return keys.keys()
    
    
    
if __name__ == "__main__":
    
    print 'Testing function : get_FilesList() : pattern: ["*.pl"] in user directory '
    files=get_FilesList('/house/homedirs/a/apratap/dev/eclipse_workspace/perl_scripts/','*.pl')
    if files:
        print 'Found %d files ' % (len(files))
    else:
        print 'No files found'
     
        
def get_full_username():
    return pwd.getpwuid(os.getuid())[4]



def merge_files(list_of_input_files,outPath='./',outFile='merged.data'):
    """
    takes a list of input files and join them into one file
    output : full path to the merged file name
    """
    
    #open the file handle for each file
    infiles_handles = ( open(f,'r') for f in list_of_input_files )
    full_out_path = outPath + '/' + outFile
    
    with open(full_out_path,'w') as fout:
        for fh in infiles_handles:
            [fout.write(line) for line in fh]
    return full_out_path



            



class NestedDict(dict):
    def __getitem__(self, key):
        if key in self: 
            return self.get(key)
        else:
            return self.setdefault(key,NestedDict())





class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value














'''
def  determineNumberOfCPUs():
    """ Number of virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError,NotImplementedError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError,ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])

        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0:
            return res
    except ImportError:
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                      stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)

        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')

        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        expr = re.compile('^cpuid@[0-9]+$')

        res = 0
        for pd in pseudoDevices:
            if expr.match(pd) != None:
                res += 1

        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0:
            return res
    except OSError:
        pass

    raise Exception('Can not determine number of CPUs on this system'
'''