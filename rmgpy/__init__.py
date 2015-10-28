#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This is the rmg module.
"""

import os
import os.path
import logging
from .version import __version__

################################################################################

class SettingsError(Exception):
    """
    An exception raised when dealing with settings.
    """
    pass

class Settings(dict):
    """
    A dictionary-like object containing global settings for RMG jobs. These
    settings are generally loaded from a file named ``rmgrc`` located at one
    of several possible places on disk. The ``rmgrc`` file used is stored in
    the `filename` attribute. This class inherits from the built-in dict class
    and adds methods for loading and resetting the settings, as well as a
    custom :meth:`__setitem__` method for processing the setting values before
    adding them to the dictionary.
    
    In general you should be working with the module-level variable
    ``settings`` in this module, which is an instance of this class.
    """
    
    def __init__(self):
        super(Settings, self).__init__()
        self.filename = None
        self.sources = dict()
        self.load()
    
    def __setitem__(self, key, value):
        if key == 'database.directory':
            value = os.path.abspath(os.path.expandvars(value))
        else:
            print('Unexpecting setting "{0}" encountered.'.format(key))
        self.sources[key] = '-'
        super(Settings, self).__setitem__(key, value)
    
    def report(self):
        """
        Returns a string saying what is set and where things came from, suitable for logging
        """
        lines = ['Global RMG Settings:']
        for key in self.keys():
            lines.append("   {0:20s} = {1:20s} ({2})".format(key, self[key], self.sources[key]))
        return '\n'.join(lines)

    def load(self, path=None):
        """
        Load settings from a file on disk. If an explicit file is not specified,
        the following locations will be searched for a settings file, and the
        first one found will be loaded:
        
        * An rmgrc file in the current working directory
        
        * An rmgrc file in the user's $HOME/.rmg directory
        
        * An rmgrc file in the same directory as this package
        
        If none of these can be found, a SettingsError is raised.
        """
        # First set all settings to their default values
        self.reset()
    
        if path:
            # The user specified an explicit file to use for the settings
            # Make sure that it exists, fail if it does not
            if not os.path.exists(path):
                raise SettingsError('Specified RMG settings file "{0}" does not exist.'.format(path))
            else:
                self.filename = path
        else:
            # The user did not specify an explicit file to use for the settings
            # Load one of the default settings files instead
            working_dir = os.path.abspath(os.path.dirname(__file__))
            if os.path.exists('rmgrc'):
                self.filename = 'rmgrc'
            elif os.path.exists(os.path.expanduser('~/.rmg/rmgrc')):
                self.filename = os.path.expanduser('~/.rmg/rmgrc')
            elif os.path.exists(os.path.join(working_dir, 'rmgrc')):
                self.filename = os.path.join(working_dir, 'rmgrc')
            else:
                return # fail silently, instead of raising the following error:
                raise SettingsError('Could not find an RMG settings file to load!')
        
        # From here on we assume that we have identified the appropriate
        # settings file to load

        with open(self.filename, 'r') as f:
            for line in f:
                # Remove any comments from the line
                index = line.find('#')
                if index != -1: line = line[:index]
                # Is there a key-value pair remaining?
                if line.find(':') != -1:
                    key, value = line.split(':')
                    key = key.strip()
                    value = value.strip()
                    self[key] = value
                    self.sources[key] = "from {0}".format(self.filename)
    
    def reset(self):
        """
        Reset all settings to their default values.
        """
        self.filename = None
        rmgpy_module_dir = os.path.abspath(os.path.dirname(__file__))
        self['database.directory'] = os.path.realpath(os.path.join(rmgpy_module_dir, '..', '..', 'RMG-database', 'input'))
        self.sources['database.directory'] = 'Default, relative to RMG-Py source code'

# The global settings object
settings = Settings()

################################################################################

def getPath():
    """
    Return the directory that this file is found in on disk.
    """
    return os.path.abspath(os.path.dirname(__file__))
