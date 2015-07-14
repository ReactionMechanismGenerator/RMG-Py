#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains the :class:`NetworkDrawer` class, used to generate a
depiction of a pressure-dependent reaction network.
"""

import numpy
import logging

from rmgpy.molecule.draw import MoleculeDrawer, createNewSurface

################################################################################

class NetworkDrawer:
    """
    This class provides functionality for drawing the potential energy surface
    for a pressure-dependent reaction network using the Cairo 2D graphics
    engine. The most common use case is simply::
    
        NetworkDrawer().draw(network, format='png', path='network.png')
    
    where ``network`` is the :class:`Network` object to draw. You can also
    pass a dict of options to the constructor to affect how the network is
    drawn.
    """
    
    def __init__(self, options=None):
        self.options = {
            'structures': True,
            'fontFamily': 'sans',
            'fontSizeNormal': 12,
            'Eunits': 'kJ/mol',
            'padding': 16,
            'wellWidth': 64,
            'wellSpacing': 64,
            'Eslope': 1.5,
            'TSwidth': 16,
            'E0offset': 0.0,
        }
        if options: self.options.update(options)
        self.clear()
    
    def clear(self):
        self.network = None
        self.left = 0.0
        self.top = 0.0
        self.right = 0.0
        self.bottom = 0.0
        self.surface = None
        self.cr = None

    def __getEnergyRange(self):
        """
        Return the minimum and maximum energy in J/mol on the potential energy
        surface.
        """
        
        E0min = self.network.isomers[0].E0
        E0max = E0min
        
        for isomer in self.network.isomers[1:]:
            E0 = isomer.E0
            if E0 < E0min: E0min = E0
            if E0 > E0max: E0max = E0
        for reactant in self.network.reactants:
            E0 = reactant.E0
            if E0 < E0min: E0min = E0
            if E0 > E0max: E0max = E0
        for product in self.network.products:
            E0 = product.E0
            if E0 < E0min: E0min = E0
            if E0 > E0max: E0max = E0
        for rxn in self.network.pathReactions:
            E0 = rxn.transitionState.conformer.E0.value_si
            if E0 < E0min: E0min = E0
            if E0 > E0max: E0max = E0
        
        return E0min, E0max

    def __useStructureForLabel(self, configuration):
        """
        Return ``True`` if the configuration should use molecular structures
        for its labels or ``False`` otherwise.
        """
        
        # Initialize with the current user option value
        useStructures = self.options['structures']
        
        # But don't use structures if one or more species in the configuration 
        # do not have structure data
        for spec in configuration.species:
            if spec.molecule is None or len(spec.molecule) == 0:
                useStructures = False
                break
            
        return useStructures
    
    def __getTextSize(self, text, padding=2, format='pdf'):
        """
        
        """
        try:
            import cairocffi as cairo
        except ImportError:
            import cairo
        
        # Use dummy surface to determine text extents
        surface = createNewSurface(format)
        cr = cairo.Context(surface)
        cr.set_font_size(self.options['fontSizeNormal'])
        extents = cr.text_extents(text)
        
        width = extents[2] + 2 * padding
        height = extents[3] + 2 * padding
        
        return [0, 0, width, height]    

    def __drawText(self, text, cr, x0, y0, padding=2):
        """
        
        """
        cr.save()
        cr.set_font_size(self.options['fontSizeNormal'])
        extents = cr.text_extents(text)
        cr.move_to(x0 - extents[0] - padding, y0 - extents[1] + padding)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.show_text(text)
        cr.restore()
        
        width = extents[2] + 2 * padding
        height = extents[3] + 2 * padding
        
        return [0, 0, width, height]    

    def __getLabelSize(self, configuration, format='pdf'):
        """
        
        """
        width = 0; height = 0; boundingRects = []
        if self.__useStructureForLabel(configuration):
            for spec in configuration.species:
                surface, cr, rect = MoleculeDrawer().draw(spec.molecule[0], format=format)
                boundingRects.append(list(rect))
        else:
            for spec in configuration.species:
                boundingRects.append(self.__getTextSize(spec.label, format=format))
        
        plusRect = self.__getTextSize('+', format=format)
        
        for rect in boundingRects:
            if width < rect[2]: width = rect[2]
            height += rect[3] + plusRect[3]
        height -= plusRect[3]
            
        return [0, 0, width, height] 

    def __drawLabel(self, configuration, cr, x0, y0, format='pdf'):
    
        boundingRect = self.__getLabelSize(configuration, format=format)
        padding = 2
        
        useStructures = self.__useStructureForLabel(configuration)
        y = y0
        for i, spec in enumerate(configuration.species):
            if i > 0:
                rect = self.__getTextSize('+', padding=padding, format=format)
                x = x0 - 0.5 * (rect[2] - boundingRect[2]) + 2 * padding
                self.__drawText('+', cr, x, y)
                y += rect[3]
            
            if useStructures:
                moleculeDrawer = MoleculeDrawer()
                cr.save()
                surf, c, rect = moleculeDrawer.draw(spec.molecule[0], format=format)
                cr.restore()
                x = x0 - 0.5 * (rect[2] - boundingRect[2])
                cr.save()
                moleculeDrawer.render(cr, offset=(x, y))
                cr.restore()
                y += rect[3]
            else:
                rect = self.__getTextSize(spec.label, padding=padding, format=format)
                x = x0 - 0.5 * (rect[2] - boundingRect[2]) + 2 * padding
                self.__drawText(spec.label, cr, x, y)
                y += rect[3]
        
        return boundingRect

    def draw(self, network, format, path=None):
        """
        Draw the potential energy surface for the given `network` as a Cairo
        surface of the given `format`. If `path` is given, the surface is
        saved to that location on disk.
        """
        try:
            import cairocffi as cairo
        except ImportError:
            try:
                import cairo
            except ImportError:
                logging.warning('Cairo not found; potential energy surface will not be drawn.')
                return

        self.network = network

        # The order of wells is as follows:
        #   - Reactant channels come first (to the left)
        #   - Isomers are in the middle
        #   - Product channels come last (to the right)
        # This is done because most people will read the PES from left to right
        wells = []
        wells.extend(network.reactants)
        wells.extend(network.isomers)
        wells.extend(network.products)

        # Generate the bounding rectangles for each configuration label
        labelRects = []
        for well in wells:
            labelRects.append(self.__getLabelSize(well, format=format))

        # Get energy range (use kJ/mol internally)
        E0min, E0max = self.__getEnergyRange()
        E0min *= 0.001; E0max *= 0.001
        
        # Drawing parameters
        padding = self.options['padding']
        wellWidth = self.options['wellWidth']
        wellSpacing = self.options['wellSpacing']
        Eslope = self.options['Eslope']
        TSwidth = self.options['TSwidth']
        
        E0_offset = self.options['E0offset'] * 0.001
        
        # Choose multiplier to convert energies to desired units (on figure only)
        Eunits = self.options['Eunits']
        try:
            Emult = {'J/mol': 1.0, 'kJ/mol': 0.001, 'cal/mol': 1.0/4.184, 'kcal/mol': 1.0/4184., 'cm^-1': 1.0/11.962}[Eunits]
        except KeyError:
            raise Exception('Invalid value "{0}" for Eunits parameter.'.format(Eunits))
            
        # Determine height required for drawing
        Eheight = self.__getTextSize('0.0', format=format)[3] + 6
        y_E0 = (E0max - 0.0) * Eslope + padding + Eheight
        height = (E0max - E0min) * Eslope + 2 * padding + Eheight + 6
        for i in range(len(wells)):
            if 0.001 * wells[i].E0 == E0min:
                height += labelRects[i][3]
                break
        
        # Determine naive position of each well (one per column)
        coordinates = numpy.zeros((len(wells), 2), numpy.float64)
        x = padding
        for i in range(len(wells)):
            well = wells[i]
            rect = labelRects[i]
            thisWellWidth = max(wellWidth, rect[2])
            E0 = 0.001 * well.E0
            y = y_E0 - E0 * Eslope
            coordinates[i] = [x + 0.5 * thisWellWidth, y]
            x += thisWellWidth + wellSpacing
        width = x + padding - wellSpacing
        
        # Determine the rectangles taken up by each well
        # We'll use this to merge columns safely so that wells don't overlap
        wellRects = []
        for i in range(len(wells)):
            l, t, w, h = labelRects[i]
            x, y = coordinates[i,:]
            if w < wellWidth: w = wellWidth
            t -= 6 + Eheight
            h += 6 + Eheight
            wellRects.append([l + x - 0.5 * w, t + y + 6, w, h])
        
        # Squish columns together from the left where possible until an isomer is encountered
        oldLeft = numpy.min(coordinates[:,0])
        Nleft = wells.index(network.isomers[0])-1
        columns = []
        for i in range(Nleft, -1, -1):
            top = wellRects[i][1]
            bottom = top + wellRects[i][3]
            for j in range(len(columns)):
                for c in columns[j]:
                    top0 = wellRects[c][1]
                    bottom0 = top + wellRects[c][3]
                    if (top >= top0 and top <= bottom0) or (top <= top0 and top0 <= bottom):
                        # Can't put it in this column
                        break
                else:
                    # Can put it in this column
                    columns[j].append(i)
                    break
            else:
                # Needs a new column
                columns.append([i])
        for column in columns:
            columnWidth = max([wellRects[c][2] for c in column])
            x = coordinates[column[0]+1,0] - 0.5 * wellRects[column[0]+1][2] - wellSpacing - 0.5 * columnWidth
            for c in column:
                delta = x - coordinates[c,0]
                wellRects[c][0] += delta
                coordinates[c,0] += delta
        newLeft = numpy.min(coordinates[:,0])
        coordinates[:,0] -= newLeft - oldLeft

        # Squish columns together from the right where possible until an isomer is encountered
        Nright = wells.index(network.isomers[-1])+1
        columns = []
        for i in range(Nright, len(wells)):
            top = wellRects[i][1]
            bottom = top + wellRects[i][3]
            for j in range(len(columns)):
                for c in columns[j]:
                    top0 = wellRects[c][1]
                    bottom0 = top0 + wellRects[c][3]
                    if (top >= top0 and top <= bottom0) or (top <= top0 and top0 <= bottom):
                        # Can't put it in this column
                        break
                else:
                    # Can put it in this column
                    columns[j].append(i)
                    break
            else:
                # Needs a new column
                columns.append([i])
        for column in columns:
            columnWidth = max([wellRects[c][2] for c in column])
            x = coordinates[column[0]-1,0] + 0.5 * wellRects[column[0]-1][2] + wellSpacing + 0.5 * columnWidth
            for c in column:
                delta = x - coordinates[c,0]
                wellRects[c][0] += delta
                coordinates[c,0] += delta

        width = max([rect[2]+rect[0] for rect in wellRects]) - min([rect[0] for rect in wellRects]) + 2 * padding
        
        # Draw to the final surface
        surface = createNewSurface(format=format, path=path, width=width, height=height)
        cr = cairo.Context(surface)
        
        # Some global settings
        cr.select_font_face("sans")
        cr.set_font_size(self.options['fontSizeNormal'])

        # Fill the background with white
        cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        cr.paint()

#        # DEBUG: Draw well bounding rectangles
#        cr.save()
#        cr.set_line_width(1.0)
#        for rect in wellRects:
#            cr.rectangle(*rect)
#            cr.set_source_rgba(0.0, 0.0, 1.0, 0.5)
#            cr.stroke()  
#        cr.restore()

        # Draw path reactions
        for rxn in network.pathReactions:
            for reac in range(len(wells)):
                if wells[reac].species == rxn.reactants:
                    break
            else:
                raise Exception
            for prod in range(len(wells)):
                if wells[prod].species == rxn.products:
                    break
            else:
                raise Exception
            E0_reac = wells[reac].E0 * 0.001 - E0_offset
            E0_prod = wells[prod].E0 * 0.001 - E0_offset
            E0_TS = rxn.transitionState.conformer.E0.value_si * 0.001 - E0_offset
            if reac < prod:
                x1, y1 = coordinates[reac,:]
                x2, y2 = coordinates[prod,:]
            else:
                x1, y1 = coordinates[prod,:]
                x2, y2 = coordinates[reac,:]
            x1 += wellSpacing / 2.0; x2 -= wellSpacing / 2.0
            if abs(E0_TS - E0_reac) > 0.1 and abs(E0_TS - E0_prod) > 0.1:
                if len(rxn.reactants) == 2:
                    if reac < prod: x0 = x1 + wellSpacing * 0.5
                    else:           x0 = x2 - wellSpacing * 0.5
                elif len(rxn.products) == 2:
                    if reac < prod: x0 = x2 - wellSpacing * 0.5
                    else:           x0 = x1 + wellSpacing * 0.5
                else:
                    x0 = 0.5 * (x1 + x2)
                y0 = y_E0 - (E0_TS + E0_offset) * Eslope
                width1 = (x0 - x1)
                width2 = (x2 - x0)
                # Draw horizontal line for TS
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.set_line_width(2.0)
                cr.move_to(x0 - TSwidth/2.0, y0)
                cr.line_to(x0+TSwidth/2.0, y0)
                cr.stroke()
                # Add background and text for energy
                E0 = "{0:.1f}".format(E0_TS * 1000. * Emult)
                extents = cr.text_extents(E0)
                x = x0 - extents[2] / 2.0; y = y0 - 6.0
                cr.rectangle(x + extents[0] - 2.0, y + extents[1] - 2.0, extents[2] + 4.0, extents[3] + 4.0)
                cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
                cr.fill()
                cr.move_to(x, y)
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.show_text(E0)
                # Draw Bezier curve connecting reactants and products through TS
                cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
                cr.set_line_width(1.0)
                cr.move_to(x1, y1)
                cr.curve_to(x1 + width1/8.0, y1,   x0 - width1/8.0 - TSwidth/2.0, y0,   x0 - TSwidth/2.0, y0)
                cr.move_to(x0 + TSwidth/2.0, y0)
                cr.curve_to(x0 + width2/8.0 + TSwidth/2.0, y0,   x2 - width2/8.0, y2,   x2, y2)
                cr.stroke()
            else:
                width = (x2 - x1)
                # Draw Bezier curve connecting reactants and products through TS
                cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
                cr.set_line_width(1.0)
                cr.move_to(x1, y1)
                cr.curve_to(x1 + width/4.0, y1,   x2 - width/4.0, y2,   x2, y2)
                cr.stroke()

        # Draw wells (after path reactions so that they are on top)
        for i, well in enumerate(wells):
            x0, y0 = coordinates[i,:]
            # Draw horizontal line for well
            cr.set_line_width(4.0)
            cr.move_to(x0 - wellWidth/2.0, y0)
            cr.line_to(x0 + wellWidth/2.0, y0)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.stroke()
            # Add background and text for energy
            E0 = well.E0 * 0.001 - E0_offset
            E0 = "{0:.1f}".format(E0 * 1000. * Emult)
            extents = cr.text_extents(E0)
            x = x0 - extents[2] / 2.0; y = y0 - 6.0
            cr.rectangle(x + extents[0] - 2.0, y + extents[1] - 2.0, extents[2] + 4.0, extents[3] + 4.0)
            cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
            cr.fill()
            cr.move_to(x, y)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.show_text(E0)
            # Draw background and text for label
            x = x0 - 0.5 * labelRects[i][2]
            y = y0 + 6
            cr.rectangle(x, y, labelRects[i][2], labelRects[i][3])
            cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
            cr.fill()
            self.__drawLabel(well, cr, x, y, format=format)
        
        # Finish Cairo drawing
        if format == 'png':
            surface.write_to_png(path)
        else:
            surface.finish()
