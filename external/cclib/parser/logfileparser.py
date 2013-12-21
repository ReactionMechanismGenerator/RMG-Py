"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 879 $"


import StringIO

try:
    import bz2 # New in Python 2.3.
except ImportError:
    bz2 = None
import fileinput
import gzip
import inspect
import logging
logging.logMultiprocessing =  0 # To avoid a problem with Avogadro
import os
import random
try:
    set # Standard type from Python 2.4+.
except NameError:
    from sets import Set as set
import sys
import types
import zipfile

import numpy

import utils
from data import ccData


def openlogfile(filename):
    """Return a file object given a filename.

    Given the filename of a log file or a gzipped, zipped, or bzipped
    log file, this function returns a regular Python file object.
    
    Given an address starting with http://, this function retrieves the url
    and returns a file object using a temporary file.

    Given a list of filenames, this function returns a FileInput object,
    which can be used for seamless iteration without concatenation.
    """

    # If there is a single string argument given.
    if type(filename) in [str, unicode]:

        extension = os.path.splitext(filename)[1]
        
        if extension == ".gz":
            fileobject = gzip.open(filename, "r")

        elif extension == ".zip":
            zip = zipfile.ZipFile(filename, "r")
            assert len(zip.namelist()) == 1, "ERROR: Zip file contains more than 1 file"
            fileobject = StringIO.StringIO(zip.read(zip.namelist()[0]))

        elif extension in ['.bz', '.bz2']:
            # Module 'bz2' is not always importable.
            assert bz2 != None, "ERROR: module bz2 cannot be imported"
            fileobject = bz2.BZ2File(filename, "r")

        else:
            fileobject = open(filename, "r")

        return fileobject
    
    elif hasattr(filename, "__iter__"):
    
        # Compression (gzip and bzip) is supported as of Python 2.5.
        if sys.version_info[0] >= 2 and sys.version_info[1] >= 5:
            fileobject = fileinput.input(filename, openhook=fileinput.hook_compressed)
        else:
            fileobject = fileinput.input(filename)
        
        return fileobject


class Logfile(object):
    """Abstract class for logfile objects.

    Subclasses defined by cclib:
        ADF, GAMESS, GAMESSUK, Gaussian, Jaguar, Molpro, ORCA
    
    """

    def __init__(self, source, progress=None,
                       loglevel=logging.INFO, logname="Log", logstream=sys.stdout,
                       fupdate=0.05, cupdate=0.002, 
                       datatype=ccData):
        """Initialise the Logfile object.

        This should be called by a ubclass in its own __init__ method.

        Inputs:
            source - a single logfile, a list of logfiles, or input stream
        """

        # Set the filename to source if it is a string or a list of filenames.
        # In the case of an input stream, set some arbitrary name and the stream.
        # Elsewise, raise an Exception.
        if isinstance(source,types.StringTypes):
            self.filename = source
            self.isstream = False
        elif isinstance(source,list) and all([isinstance(s,types.StringTypes) for s in source]):
            self.filename = source
            self.isstream = False
        elif hasattr(source, "read"):
            self.filename = "stream %s" %str(type(source))
            self.isstream = True
            self.stream = source
        else:
            raise ValueError

        # Progress indicator.
        self.progress = progress
        self.fupdate = fupdate
        self.cupdate = cupdate

        # Set up the logger.
        # Note that calling logging.getLogger() with one name always returns the same instance.
        # Presently in cclib, all parser instances of the same class use the same logger,
        #   which means that care needs to be taken not to duplicate handlers.
        self.loglevel = loglevel
        self.logname  = logname
        self.logger = logging.getLogger('%s %s' % (self.logname,self.filename))
        self.logger.setLevel(self.loglevel)
        if len(self.logger.handlers) == 0:
                handler = logging.StreamHandler(logstream)
                handler.setFormatter(logging.Formatter("[%(name)s %(levelname)s] %(message)s"))
                self.logger.addHandler(handler)

        # Periodic table of elements.
        self.table = utils.PeriodicTable()

        # This is the class that will be used in the data object returned by parse(),
        #   and should normally be ccData or a subclass.
        self.datatype = datatype

    def __setattr__(self, name, value):

        # Send info to logger if the attribute is in the list self._attrlist.
        if name in getattr(self, "_attrlist", {}) and hasattr(self, "logger"):
                    
            # Call logger.info() only if the attribute is new.
            if not hasattr(self, name):
                if type(value) in [numpy.ndarray, list]:
                    self.logger.info("Creating attribute %s[]" %name)
                else:
                    self.logger.info("Creating attribute %s: %s" %(name, str(value)))

        # Set the attribute.
        object.__setattr__(self, name, value)

    def parse(self, fupdate=None, cupdate=None):
        """Parse the logfile, using the assumed extract method of the child."""

        # Check that the sub-class has an extract attribute,
        #  that is callable with the proper number of arguemnts.
        if not hasattr(self, "extract"):
            raise AttributeError, "Class %s has no extract() method." %self.__class__.__name__
            return -1
        if not callable(self.extract):
            raise AttributeError, "Method %s._extract not callable." %self.__class__.__name__
            return -1
        if len(inspect.getargspec(self.extract)[0]) != 3:
            raise AttributeError, "Method %s._extract takes wrong number of arguments." %self.__class__.__name__
            return -1

        # Save the current list of attributes to keep after parsing.
        # The dict of self should be the same after parsing.
        _nodelete = list(set(self.__dict__.keys()))

        # Initiate the FileInput object for the input files.
        # Remember that self.filename can be a list of files.
        if not self.isstream:
            inputfile = openlogfile(self.filename)
        else:
            inputfile = self.stream

        # Intialize self.progress.
        if self.progress:
            inputfile.seek(0,2)
            nstep = inputfile.tell()
            inputfile.seek(0)
            self.progress.initialize(nstep)
            self.progress.step = 0
            if fupdate:
                self.fupdate = fupdate
            if cupdate:
                self.cupdate = cupdate

        # Initialize the ccData object that will be returned.
        # This is normally ccData, but can be changed by passing
        #   the datatype argument to __init__().
        data = self.datatype()
        
        # Copy the attribute list, so that the parser knows what to expect,
        #   specifically in __setattr__().
        # The class self.datatype (normally ccData) must have this attribute.
        self._attrlist = data._attrlist
        
        # Maybe the sub-class has something to do before parsing.
        if hasattr(self, "before_parsing"):
            self.before_parsing()

        # Loop over lines in the file object and call extract().
        # This is where the actual parsing is done.
        for line in inputfile:

            self.updateprogress(inputfile, "Unsupported information", cupdate)

            # This call should check if the line begins a section of extracted data.
            # If it does, it parses some lines and sets the relevant attributes (to self).
            # Any attributes can be freely set and used across calls, however only those
            #   in data._attrlist will be moved to final data object that is returned.
            self.extract(inputfile, line)

        # Close input file object.
        if not self.isstream:
            inputfile.close()

        # Maybe the sub-class has something to do after parsing.
        if hasattr(self, "after_parsing"):
            self.after_parsing()

        # If atomcoords were not parsed, but some input coordinates were ("inputcoords").
        # This is originally from the Gaussian parser, a regression fix.
        if not hasattr(self, "atomcoords") and hasattr(self, "inputcoords"):
            self.atomcoords = numpy.array(self.inputcoords, 'd')

        # Set nmo if not set already - to nbasis.
        if not hasattr(self, "nmo") and hasattr(self, "nbasis"):
            self.nmo = self.nbasis

        # Creating deafult coreelectrons array.
        if not hasattr(self, "coreelectrons") and hasattr(self, "natom"):
            self.coreelectrons = numpy.zeros(self.natom, "i")

        # Move all cclib attributes to the ccData object.
        # To be moved, an attribute must be in data._attrlist.
        for attr in data._attrlist:
            if hasattr(self, attr):
                setattr(data, attr, getattr(self, attr))
                
        # Now make sure that the cclib attributes in the data object
        #   are all the correct type (including arrays and lists of arrays).
        data.arrayify()

        # Delete all temporary attributes (including cclib attributes).
        # All attributes should have been moved to a data object,
        #   which will be returned.
        for attr in self.__dict__.keys():
            if not attr in _nodelete:
                self.__delattr__(attr)

        # Update self.progress as done.
        if self.progress:
            self.progress.update(nstep, "Done")

        # Return the ccData object that was generated.
        return data

    def updateprogress(self, inputfile, msg, xupdate=0.05):
        """Update progress."""

        if self.progress and random.random() < xupdate:
            newstep = inputfile.tell()
            if newstep != self.progress.step:
                self.progress.update(newstep, msg)
                self.progress.step = newstep

    def normalisesym(self,symlabel):
        """Standardise the symmetry labels between parsers.

        This method should be overwritten by individual parsers, and should
        contain appropriate doctests. If is not overwritten, this is detected
        as an error by unit tests.
        """
        return "ERROR: This should be overwritten by this subclass"

    def float(self,number):
        """Convert a string to a float avoiding the problem with Ds.

        >>> t = Logfile("dummyfile")
        >>> t.float("123.2323E+02")
        12323.23
        >>> t.float("123.2323D+02")
        12323.23
        """
        number = number.replace("D","E")
        return float(number)

if __name__=="__main__":
    import doctest
    doctest.testmod()
