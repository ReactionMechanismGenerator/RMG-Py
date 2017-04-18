import os
import os.path
import sys
import logging
from functools import wraps
import cProfile
import pstats
import subprocess

class Tee:
    """A simple tee to create a stream which prints to many streams.
    
    This is used to report the profiling statistics to both the log file
    and the standard output.
    """
    def __init__(self, *fileobjects):
        self.fileobjects=fileobjects
    def write(self, string):
        for fileobject in self.fileobjects:
            fileobject.write(string)

def processProfileStats(stats_file, log_file):
    
    out_stream = Tee(sys.stdout,open(log_file,'a')) # print to screen AND append to RMG.log
    print >>out_stream, "="*80
    print >>out_stream, "Profiling Data".center(80)
    print >>out_stream, "="*80
    stats = pstats.Stats(stats_file,stream=out_stream)
    stats.strip_dirs()
    print >>out_stream, "Sorted by internal time"
    stats.sort_stats('time')
    stats.print_stats(25)
    stats.print_callers(25)
    print >>out_stream, "Sorted by cumulative time"
    stats.sort_stats('cumulative')
    stats.print_stats(25)
    stats.print_callers(25)
    stats.print_callees(25)

def makeProfileGraph(stats_file):
    """
    Uses gprof2dot to create a graphviz dot file of the profiling information.
    
    This requires the gprof2dot package available via `pip install gprof2dot`.
    Render the result using the program 'dot' via a command like
    `dot -Tps2 input.dot -o output.ps2`.
    
    Rendering the ps2 file to pdf requires an external pdf converter
    `ps2pdf output.ps2` which produces a `output.ps2.pdf` file.
    """
    try:
        from gprof2dot import PstatsParser, DotWriter, SAMPLES, themes
    except ImportError:
        logging.warning('Trouble importing from package gprof2dot. Unable to create a graph of the profile statistics.')
        logging.warning('Try getting the latest version with something like `pip install --upgrade gprof2dot`.')
        return
    
    #create an Options class to mimic optparser output as much as possible:
    class Options:
        pass
    
    options = Options()
    options.node_thres = 0.8
    options.edge_thres = 0.1
    options.strip = False
    options.show_samples = False
    options.root = ""
    options.leaf = ""
    options.wrap = True
    
    theme = themes['color'] # bw color gray pink
    theme.fontname = "ArialMT" # default "Arial" leads to PostScript warnings in dot (on Mac OS)
    parser = PstatsParser(stats_file)
    profile = parser.parse()
    
    dot_file = stats_file + '.dot'
    output = open(dot_file,'wt')
    dot = DotWriter(output)
    dot.strip = options.strip
    dot.wrap = options.wrap
    
    if options.show_samples:
        dot.show_function_events.append(SAMPLES)
    
    profile = profile
    profile.prune(options.node_thres/100.0, options.edge_thres/100.0)

    if options.root:
        rootId = profile.getFunctionId(options.root)
        if not rootId:
            sys.stderr.write('root node ' + options.root + ' not found (might already be pruned : try -e0 -n0 flags)\n')
            sys.exit(1)
        profile.prune_root(rootId)
    if options.leaf:
        leafId = profile.getFunctionId(options.leaf)
        if not leafId:
            sys.stderr.write('leaf node ' + options.leaf + ' not found (maybe already pruned : try -e0 -n0 flags)\n')
            sys.exit(1)
        profile.prune_leaf(leafId)

    dot.graph(profile, theme)

    output.close()
    
    try:
        subprocess.check_call(['dot', '-Tps2', dot_file, '-o', '{0}.ps2'.format(dot_file)])
    except subprocess.CalledProcessError:
        logging.error("Error returned by 'dot' when generating graph of the profile statistics.")
        logging.info("To try it yourself:\n     dot -Tps2 {0} -o {0}.ps2".format(dot_file))
    except OSError:
        logging.error("Couldn't run 'dot' to create graph of profile statistics. Check graphviz is installed properly and on your path.")
        logging.info("Once you've got it, try:\n     dot -Tps2 {0} -o {0}.ps2".format(dot_file))
    
    try:
        subprocess.check_call(['ps2pdf', '{0}.ps2'.format(dot_file), '{0}.pdf'.format(dot_file)])
    except OSError:
        logging.error("Couldn't run 'ps2pdf' to create pdf graph of profile statistics. Check that ps2pdf converter is installed.")
        logging.info("Once you've got it, try:\n     pd2pdf {0}.ps2 {0}.pdf".format(dot_file))    
    else:
        logging.info("Graph of profile statistics saved to: \n {0}.pdf".format(dot_file))

class profiler(object):
   "Decorator that keeps track of the number of times a function is called."

   __instances = {}

   def __init__(self, f):
      self.__f = f
      self.__numcalls = 0
      profiler.__instances[f] = self

   def __call__(self, *args, **kwargs):
      self.__numcalls += 1
      with self:
        return self.profile(*args, **kwargs)

   def __enter__(self):
        return self

   def __exit__(self, exc_type, exc_value, traceback):
        module = sys.modules[self.__f.__module__]
        dirname = os.path.join(os.path.dirname(module.__file__))

        stats_file = os.path.join(dirname, self.__f.__name__+'.profile')

        stats_i = os.path.join(
            dirname, ''.join([self.__f.__name__, str(self.count()), '.profile'])
                )

        if self.count() == 1:
            stats = pstats.Stats(stats_i)
        else:
            stats = pstats.Stats(stats_file)
            stats.add(stats_i)
        
        stats.dump_stats(stats_file)
        os.unlink(stats_i)

        makeProfileGraph(stats_file)

   def count(self):
      "Return the number of times the function f was called."
      return profiler.__instances[self.__f].__numcalls

   @staticmethod
   def counts():
      "Return a dict of {function: # of calls} for all registered functions."
      return dict([(f.__name__, profiler.__instances[f].__numcalls) for f in profiler.__instances])

   def profile(self, *args, **kwargs):
        module = sys.modules[self.__f.__module__]
        dirname = os.path.join(os.path.dirname(module.__file__))

        prof = cProfile.Profile()
        retval = prof.runcall(self.__f, *args, **kwargs)

        stats_file = os.path.join(
            dirname, ''.join([self.__f.__name__, str(self.count()), '.profile'])
                )

        prof.dump_stats(stats_file)
        return retval