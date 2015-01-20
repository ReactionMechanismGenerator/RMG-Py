# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

import logging
import sys

"""Abstract based method class."""

class Method(object):
    """Abstract class for method classes.

    Subclasses defined by cclib:
        CDA, CSPA, Density, FragmentAnalysis, LPA, MBO, MPA, Nuclear, OPA, Population, Volume

    All the modules containing methods should be importable:
    >>> import cda, cspa, density, fragments, lpa, mbo, mpa, nuclear, opa, population, volume
    """

    def __init__(self, data, progress=None, loglevel=logging.INFO, logname="Log"):
        """Initialise the Logfile object.

        This constructor is typically called by the constructor of a subclass.
        """

        self.data = data
        self.progress = progress
        self.loglevel = loglevel
        self.logname = logname

        self.logger = logging.getLogger('%s %s' % (self.logname, self.data))
        self.logger.setLevel(self.loglevel)
        self.logformat = "[%(name)s %(levelname)s] %(message)s"
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(self.logformat))
        self.logger.addHandler(handler)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=False)
