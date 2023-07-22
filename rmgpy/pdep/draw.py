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

"""
This module contains the :class:`NetworkDrawer` class, used to generate a
depiction of a pressure-dependent reaction network.
"""

import logging

import numpy as np

from rmgpy.molecule.draw import MoleculeDrawer, create_new_surface


################################################################################


class NetworkDrawer(object):
    """
    This class provides functionality for drawing the potential energy surface
    for a pressure-dependent reaction network using the Cairo 2D graphics
    engine. The most common use case is simply::

        NetworkDrawer().draw(network, file_format='png', path='network.png')

    where ``network`` is the :class:`Network` object to draw. You can also
    pass a dict of options to the constructor to affect how the network is
    drawn.
    """

    def __init__(self, options=None):
        self.options = {
            "structures": True,
            "fontFamily": "sans",
            "fontSizeNormal": 12,
            "Eunits": "kJ/mol",
            "padding": 16,
            "wellWidth": 64,
            "wellSpacing": 64,
            "Eslope": 1.5,
            "TSwidth": 16,
            "E0offset": 0.0,
        }
        if options:
            self.options.update(options)

        self.network = None
        self.left = 0.0
        self.top = 0.0
        self.right = 0.0
        self.bottom = 0.0
        self.surface = None
        self.cr = None

    def clear(self):
        self.network = None
        self.left = 0.0
        self.top = 0.0
        self.right = 0.0
        self.bottom = 0.0
        self.surface = None
        self.cr = None

    def _get_energy_range(self):
        """
        Return the minimum and maximum energy in J/mol on the potential energy
        surface.
        """

        e0_min = self.network.isomers[0].E0
        e0_max = e0_min

        for isomer in self.network.isomers[1:]:
            E0 = isomer.E0
            if E0 < e0_min:
                e0_min = E0
            if E0 > e0_max:
                e0_max = E0
        for reactant in self.network.reactants:
            E0 = reactant.E0
            if E0 < e0_min:
                e0_min = E0
            if E0 > e0_max:
                e0_max = E0
        for product in self.network.products:
            E0 = product.E0
            if E0 < e0_min:
                e0_min = E0
            if E0 > e0_max:
                e0_max = E0
        for rxn in self.network.path_reactions:
            E0 = rxn.transition_state.conformer.E0.value_si
            if E0 < e0_min:
                e0_min = E0
            if E0 > e0_max:
                e0_max = E0

        return e0_min, e0_max

    def _use_structure_for_label(self, configuration):
        """
        Return ``True`` if the configuration should use molecular structures
        for its labels or ``False`` otherwise.
        """

        # Initialize with the current user option value
        use_structures = self.options["structures"]

        # But don't use structures if one or more species in the configuration
        # do not have structure data
        for spec in configuration.species:
            if spec.molecule is None or len(spec.molecule) == 0:
                use_structures = False
                break

        return use_structures

    def _get_text_size(self, text, padding=2, file_format="pdf"):
        """ """
        try:
            import cairocffi as cairo
        except ImportError:
            import cairo

        # Use dummy surface to determine text extents
        surface = create_new_surface(file_format)
        cr = cairo.Context(surface)
        cr.set_font_size(self.options["fontSizeNormal"])
        extents = cr.text_extents(text)

        width = extents[2] + 2 * padding
        height = extents[3] + 2 * padding

        return [0, 0, width, height]

    def _draw_text(self, text, cr, x0, y0, padding=2):
        """ """
        cr.save()
        cr.set_font_size(self.options["fontSizeNormal"])
        extents = cr.text_extents(text)
        cr.move_to(x0 - extents[0] - padding, y0 - extents[1] + padding)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.show_text(text)
        cr.restore()

        width = extents[2] + 2 * padding
        height = extents[3] + 2 * padding

        return [0, 0, width, height]

    def _get_label_size(self, configuration, file_format="pdf"):
        """ """
        width = 0
        height = 0
        bounding_rects = []
        if self._use_structure_for_label(configuration):
            for spec in configuration.species:
                surface, cr, rect = MoleculeDrawer().draw(spec.molecule[0], file_format=file_format)
                bounding_rects.append(list(rect))
        else:
            for spec in configuration.species:
                bounding_rects.append(self._get_text_size(spec.label, file_format=file_format))

        plus_rect = self._get_text_size("+", file_format=file_format)

        for rect in bounding_rects:
            if width < rect[2]:
                width = rect[2]
            height += rect[3] + plus_rect[3]
        height -= plus_rect[3]

        return [0, 0, width, height]

    def _draw_label(self, configuration, cr, x0, y0, file_format="pdf"):
        bounding_rect = self._get_label_size(configuration, file_format=file_format)
        padding = 2

        use_structures = self._use_structure_for_label(configuration)
        y = y0
        for i, spec in enumerate(configuration.species):
            if i > 0:
                rect = self._get_text_size("+", padding=padding, file_format=file_format)
                x = x0 - 0.5 * (rect[2] - bounding_rect[2]) + 2 * padding
                self._draw_text("+", cr, x, y)
                y += rect[3]

            if use_structures:
                molecule_drawer = MoleculeDrawer()
                cr.save()
                surf, c, rect = molecule_drawer.draw(spec.molecule[0], file_format=file_format)
                cr.restore()
                x = x0 - 0.5 * (rect[2] - bounding_rect[2])
                cr.save()
                molecule_drawer.render(cr, offset=(x, y))
                cr.restore()
                y += rect[3]
            else:
                rect = self._get_text_size(spec.label, padding=padding, file_format=file_format)
                x = x0 - 0.5 * (rect[2] - bounding_rect[2]) + 2 * padding
                self._draw_text(spec.label, cr, x, y)
                y += rect[3]

        return bounding_rect

    def draw(self, network, file_format, path=None):
        """
        Draw the potential energy surface for the given `network` as a Cairo
        surface of the given `file_format`. If `path` is given, the surface is
        saved to that location on disk.
        """
        try:
            import cairocffi as cairo
        except ImportError:
            try:
                import cairo
            except ImportError:
                logging.warning("Cairo not found; potential energy surface will not be drawn.")
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
        label_rects = []
        for well in wells:
            label_rects.append(self._get_label_size(well, file_format=file_format))

        # Get energy range (use kJ/mol internally)
        e0_min, e0_max = self._get_energy_range()
        e0_min *= 0.001
        e0_max *= 0.001

        # Drawing parameters
        padding = self.options["padding"]
        well_width = self.options["wellWidth"]
        well_spacing = self.options["wellSpacing"]
        e_slope = self.options["Eslope"]
        ts_width = self.options["TSwidth"]

        e0_offset = self.options["E0offset"] * 0.001

        # Choose multiplier to convert energies to desired units (on figure only)
        e_units = self.options["Eunits"]
        try:
            e_mult = {
                "J/mol": 1.0,
                "kJ/mol": 0.001,
                "cal/mol": 1.0 / 4.184,
                "kcal/mol": 1.0 / 4184.0,
                "cm^-1": 1.0 / 11.962,
            }[e_units]
        except KeyError:
            raise Exception('Invalid value "{0}" for Eunits parameter.'.format(e_units))

        # Determine height required for drawing
        e_height = self._get_text_size("0.0", file_format=file_format)[3] + 6
        y_e0 = (e0_max - 0.0) * e_slope + padding + e_height
        height = (e0_max - e0_min) * e_slope + 2 * padding + e_height + 6
        for i in range(len(wells)):
            if 0.001 * wells[i].E0 == e0_min:
                height += label_rects[i][3]
                break

        # Determine naive position of each well (one per column)
        coordinates = np.zeros((len(wells), 2), float)
        x = padding
        for i in range(len(wells)):
            well = wells[i]
            rect = label_rects[i]
            this_well_width = max(well_width, rect[2])
            E0 = 0.001 * well.E0
            y = y_e0 - E0 * e_slope
            coordinates[i] = [x + 0.5 * this_well_width, y]
            x += this_well_width + well_spacing
        width = x + padding - well_spacing

        # Determine the rectangles taken up by each well
        # We'll use this to merge columns safely so that wells don't overlap
        well_rects = []
        for i in range(len(wells)):
            l, t, w, h = label_rects[i]
            x, y = coordinates[i, :]
            if w < well_width:
                w = well_width
            t -= 6 + e_height
            h += 6 + e_height
            well_rects.append([l + x - 0.5 * w, t + y + 6, w, h])

        # Squish columns together from the left where possible until an isomer is encountered
        old_left = np.min(coordinates[:, 0])
        n_left = wells.index(network.isomers[0]) - 1
        columns = []
        for i in range(n_left, -1, -1):
            top = well_rects[i][1]
            bottom = top + well_rects[i][3]
            for j in range(len(columns)):
                for c in columns[j]:
                    top0 = well_rects[c][1]
                    bottom0 = top + well_rects[c][3]
                    if (top0 <= top <= bottom0) or (top <= top0 <= bottom):
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
            column_width = max([well_rects[c][2] for c in column])
            x = coordinates[column[0] + 1, 0] - 0.5 * well_rects[column[0] + 1][2] - well_spacing - 0.5 * column_width
            for c in column:
                delta = x - coordinates[c, 0]
                well_rects[c][0] += delta
                coordinates[c, 0] += delta
        new_left = np.min(coordinates[:, 0])
        coordinates[:, 0] -= new_left - old_left

        # Squish columns together from the right where possible until an isomer is encountered
        n_right = wells.index(network.isomers[-1]) + 1
        columns = []
        for i in range(n_right, len(wells)):
            top = well_rects[i][1]
            bottom = top + well_rects[i][3]
            for j in range(len(columns)):
                for c in columns[j]:
                    top0 = well_rects[c][1]
                    bottom0 = top0 + well_rects[c][3]
                    if (top0 <= top <= bottom0) or (top <= top0 <= bottom):
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
            column_width = max([well_rects[c][2] for c in column])
            x = coordinates[column[0] - 1, 0] + 0.5 * well_rects[column[0] - 1][2] + well_spacing + 0.5 * column_width
            for c in column:
                delta = x - coordinates[c, 0]
                well_rects[c][0] += delta
                coordinates[c, 0] += delta

        width = max([rect[2] + rect[0] for rect in well_rects]) - min([rect[0] for rect in well_rects]) + 2 * padding

        # Draw to the final surface
        surface = create_new_surface(file_format=file_format, target=path, width=width, height=height)
        cr = cairo.Context(surface)

        # Some global settings
        cr.select_font_face("sans")
        cr.set_font_size(self.options["fontSizeNormal"])

        # Fill the background with white
        cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        cr.paint()
        self._draw_text("E0 ({0})".format(e_units), cr, 15, 10, padding=2)  # write units

        # # DEBUG: Draw well bounding rectangles
        # cr.save()
        # cr.set_line_width(1.0)
        # for rect in wellRects:
        #    cr.rectangle(*rect)
        #    cr.set_source_rgba(0.0, 0.0, 1.0, 0.5)
        #    cr.stroke()
        # cr.restore()

        # Draw path reactions
        for rxn in network.path_reactions:
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
            e0_reac = wells[reac].E0 * 0.001 - e0_offset
            e0_prod = wells[prod].E0 * 0.001 - e0_offset
            e0_ts = rxn.transition_state.conformer.E0.value_si * 0.001 - e0_offset
            if reac < prod:
                x1, y1 = coordinates[reac, :]
                x2, y2 = coordinates[prod, :]
            else:
                x1, y1 = coordinates[prod, :]
                x2, y2 = coordinates[reac, :]
            x1 += well_spacing / 2.0
            x2 -= well_spacing / 2.0
            if abs(e0_ts - e0_reac) > 0.1 and abs(e0_ts - e0_prod) > 0.1:
                if len(rxn.reactants) == 2:
                    if reac < prod:
                        x0 = x1 + well_spacing * 0.5
                    else:
                        x0 = x2 - well_spacing * 0.5
                elif len(rxn.products) == 2:
                    if reac < prod:
                        x0 = x2 - well_spacing * 0.5
                    else:
                        x0 = x1 + well_spacing * 0.5
                else:
                    x0 = 0.5 * (x1 + x2)
                y0 = y_e0 - (e0_ts + e0_offset) * e_slope
                width1 = x0 - x1
                width2 = x2 - x0
                # Draw horizontal line for TS
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.set_line_width(2.0)
                cr.move_to(x0 - ts_width / 2.0, y0)
                cr.line_to(x0 + ts_width / 2.0, y0)
                cr.stroke()
                # Add background and text for energy
                E0 = "{0:.1f}".format(e0_ts * 1000.0 * e_mult)
                extents = cr.text_extents(E0)
                x = x0 - extents[2] / 2.0
                y = y0 - 6.0
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
                cr.curve_to(x1 + width1 / 8.0, y1, x0 - width1 / 8.0 - ts_width / 2.0, y0, x0 - ts_width / 2.0, y0)
                cr.move_to(x0 + ts_width / 2.0, y0)
                cr.curve_to(x0 + width2 / 8.0 + ts_width / 2.0, y0, x2 - width2 / 8.0, y2, x2, y2)
                cr.stroke()
            else:
                width = x2 - x1
                # Draw Bezier curve connecting reactants and products through TS
                cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
                cr.set_line_width(1.0)
                cr.move_to(x1, y1)
                cr.curve_to(x1 + width / 4.0, y1, x2 - width / 4.0, y2, x2, y2)
                cr.stroke()

        # Draw wells (after path reactions so that they are on top)
        for i, well in enumerate(wells):
            x0, y0 = coordinates[i, :]
            # Draw horizontal line for well
            cr.set_line_width(4.0)
            cr.move_to(x0 - well_width / 2.0, y0)
            cr.line_to(x0 + well_width / 2.0, y0)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.stroke()
            # Add background and text for energy
            E0 = well.E0 * 0.001 - e0_offset
            E0 = "{0:.1f}".format(E0 * 1000.0 * e_mult)
            extents = cr.text_extents(E0)
            x = x0 - extents[2] / 2.0
            y = y0 - 6.0
            cr.rectangle(x + extents[0] - 2.0, y + extents[1] - 2.0, extents[2] + 4.0, extents[3] + 4.0)
            cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
            cr.fill()
            cr.move_to(x, y)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.show_text(E0)
            # Draw background and text for label
            x = x0 - 0.5 * label_rects[i][2]
            y = y0 + 6
            cr.rectangle(x, y, label_rects[i][2], label_rects[i][3])
            cr.set_source_rgba(1.0, 1.0, 1.0, 0.75)
            cr.fill()
            self._draw_label(well, cr, x, y, file_format=file_format)

        # Finish Cairo drawing
        if file_format == "png":
            surface.write_to_png(path)
        else:
            surface.finish()
