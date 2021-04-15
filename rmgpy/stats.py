#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import logging
import os.path

import matplotlib.pyplot as plt
try:
    import xlwt
except ImportError:
    logging.warning('Optional package dependency "xlwt" not loaded. Some output features will not work.')
    xlwt = None

from rmgpy.util import make_output_subdirectory


class ExecutionStatsWriter(object):
    """
    This class listens to a RMG subject
    and writes an excel file with the memory footprint
    requirements through the course of an RMG simulation,
    
    It also generates a number of images with information on the core/edge
    species/reaction evolutions through the course of an RMG simulation.

    Files are written to the 'plot' subfolder.


    A new instance of the class can be appended to a subject as follows:
    
    rmg = ...
    listener = ExecutionStatsWriter()
    rmg.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    rmg.detach(listener)
    
    """

    def __init__(self, outputDirectory):
        super(ExecutionStatsWriter, self).__init__()
        make_output_subdirectory(outputDirectory, 'plot')

        # RMG execution statistics
        self.coreSpeciesCount = []
        self.coreReactionCount = []
        self.edgeSpeciesCount = []
        self.edgeReactionCount = []
        self.memoryUse = []

    def update(self, rmg):
        self.update_execution(rmg)

    def update_execution(self, rmg):

        # Update RMG execution statistics
        logging.info('Updating RMG execution statistics...')
        core_spec, core_reac, edge_spec, edge_reac = rmg.reaction_model.get_model_size()
        self.coreSpeciesCount.append(core_spec)
        self.coreReactionCount.append(core_reac)
        self.edgeSpeciesCount.append(edge_spec)
        self.edgeReactionCount.append(edge_reac)
        elapsed = rmg.exec_time[-1]
        seconds = elapsed % 60
        minutes = (elapsed - seconds) % 3600 / 60
        hours = (elapsed - seconds - minutes * 60) % (3600 * 24) / 3600
        days = (elapsed - seconds - minutes * 60 - hours * 3600) / (3600 * 24)
        logging.info('    Execution time (DD:HH:MM:SS): '
                     '{0:02}:{1:02}:{2:02}:{3:02}'.format(int(days), int(hours), int(minutes), int(seconds)))
        try:
            import psutil
            process = psutil.Process(os.getpid())
            memory_info = process.memory_info()
            self.memoryUse.append(memory_info.rss / 1.0e6)
            logging.info('    Memory used: %.2f MB' % (self.memoryUse[-1]))
        except:
            logging.info('    Memory used: memory usage was unable to be logged')
            self.memoryUse.append(0.0)

        self.save_execution_statistics(rmg)
        if rmg.generate_plots:
            self.generate_execution_plots(rmg)

        logging.info('')

    def save_execution_statistics(self, rmg):
        """
        Save the statistics of the RMG job to an Excel spreadsheet for easy viewing
        after the run is complete. The statistics are saved to the file
        `statistics.xls` in the output directory. The ``xlwt`` package is used to
        create the spreadsheet file; if this package is not installed, no file is
        saved.
        """

        # Attempt to import the xlwt package; return if not installed
        if xlwt is None:
            logging.warning('Package xlwt not loaded. Unable to save execution statistics.')
            return

        # Create workbook and sheet for statistics to be places
        workbook = xlwt.Workbook()
        sheet = workbook.add_sheet('Statistics')

        # First column is execution time
        sheet.write(0, 0, 'Execution time (s)')
        for i, etime in enumerate(rmg.exec_time):
            sheet.write(i + 1, 0, etime)

        # Second column is number of core species
        sheet.write(0, 1, 'Core species')
        for i, count in enumerate(self.coreSpeciesCount):
            sheet.write(i + 1, 1, count)

        # Third column is number of core reactions
        sheet.write(0, 2, 'Core reactions')
        for i, count in enumerate(self.coreReactionCount):
            sheet.write(i + 1, 2, count)

        # Fourth column is number of edge species
        sheet.write(0, 3, 'Edge species')
        for i, count in enumerate(self.edgeSpeciesCount):
            sheet.write(i + 1, 3, count)

        # Fifth column is number of edge reactions
        sheet.write(0, 4, 'Edge reactions')
        for i, count in enumerate(self.edgeReactionCount):
            sheet.write(i + 1, 4, count)

        # Sixth column is memory used
        sheet.write(0, 5, 'Memory used (MB)')
        for i, memory in enumerate(self.memoryUse):
            sheet.write(i + 1, 5, memory)

        # Save workbook to file
        fstr = os.path.join(rmg.output_directory, 'statistics.xls')
        workbook.save(fstr)

    def generate_execution_plots(self, rmg):
        """
        Generate a number of plots describing the statistics of the RMG job,
        including the reaction model core and edge size and memory use versus
        execution time. These will be placed in the output directory in the plot/
        folder.
        """

        logging.info('Generating plots of execution statistics...')

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.semilogx(rmg.exec_time, self.coreSpeciesCount, 'o-b')
        ax1.set_xlabel('Execution time (s)')
        ax1.set_ylabel('Number of core species')
        ax2 = ax1.twinx()
        ax2.semilogx(rmg.exec_time, self.coreReactionCount, 'o-r')
        ax2.set_ylabel('Number of core reactions')
        plt.savefig(os.path.join(rmg.output_directory, 'plot/coreSize.svg'))
        plt.close()

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        if any(self.edgeSpeciesCount):
            ax1.loglog(rmg.exec_time, self.edgeSpeciesCount, 'o-b')
        else:
            ax1.semilogx(rmg.exec_time, self.edgeSpeciesCount, 'o-b')
        ax1.set_xlabel('Execution time (s)')
        ax1.set_ylabel('Number of edge species')
        ax2 = ax1.twinx()
        if any(self.edgeReactionCount):
            ax2.loglog(rmg.exec_time, self.edgeReactionCount, 'o-r')
        else:
            ax2.semilogx(rmg.exec_time, self.edgeReactionCount, 'o-r')
        ax2.set_ylabel('Number of edge reactions')
        plt.savefig(os.path.join(rmg.output_directory, 'plot/edgeSize.svg'))
        plt.close()

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.semilogx(rmg.exec_time, self.memoryUse, 'o-k')
        ax1.set_xlabel('Execution time (s)')
        ax1.set_ylabel('Memory (MB)')
        ax1.legend(['RAM'], loc=2)
        plt.savefig(os.path.join(rmg.output_directory, 'plot/memoryUse.svg'))
        plt.close()
