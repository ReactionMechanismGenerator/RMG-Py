#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
Contains classes and functions for working with the bibliographic information.
Currently there are three such classes:

* :class:`Article` - For articles in a journal, magazine, or other periodical

* :class:`Book` - For complete books

* :class:`Thesis` - For a graduate thesis

The above are all derived from the base :class:`Reference` class, which can
also be used if the reference does not fit into any of the above categories.
"""

################################################################################

class Reference:
    """
    A base class for representing bibliographic information. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `authors`           A list of the authors of the reference
    `title`             The title of the reference
    `year`              The year the reference was published (as a string)
    `url`               A DOI or other permalink to the reference
    =================== ========================================================

    """

    def __init__(self, authors=None, title='', year='', url=''):
        self.authors = authors or []
        self.title = title
        self.year = year
        self.url = url
    
    def getAuthorString(self):
        """
        Return a pretty, reStructuredText-formatted string of the authors.
        """
        if self.authors is not None and len(self.authors) > 0:
            if len(self.authors) == 1:
                return '%s.' % (self.authors[0])
            elif len(self.authors) == 2:
                return '%s and %s.' % (self.authors[0])
            elif self.authors[-1] == 'et al':
                return '%s et al.' % (', '.join(self.authors[:-1]))
            else:
                return '%s, and %s.' % (', '.join(self.authors[:-1]), self.authors[-1])
        else:
            return ''

################################################################################

class Article(Reference):
    """
    A class for representing an article in a journal, magazine, or other
    periodical. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `authors`           A list of the authors of the reference
    `title`             The title of the reference
    `journal`           The abbreviated name of the journal
    `volume`            The volume that the article appears in (as a string)
    `number`            The number that the article appears in (as a string)
    `pages`             The range of pages of the article (as a string)
    `year`              The year the reference was published (as a string)
    `url`               A DOI or other permalink to the reference
    =================== ========================================================

    """

    def __init__(self, authors=None, title='', journal='', volume='', number='', pages='', year='', url=''):
        Reference.__init__(self, authors=authors, title=title, year=year, url=url)
        self.journal = journal
        self.volume = volume
        self.number = number
        self.pages = pages

    def __str__(self):
        """
        Return a string representation of the reference in reStructuredText
        format.
        """
        string = self.getAuthorString()
        if self.title != '':
            string += ' "%s."' % (self.title)
        if self.journal != '':
            string += ' *%s*' % (self.journal)
        if self.volume != '':
            string += ' **%s**' % (self.volume)
        if self.number != '':
            string += ' (%s)' % (self.number)
        if self.pages != '':
            string += ', p. %s' % (self.pages)
        if self.year != '':
            string += ' (%s)' % (self.year)
        if string[-1] != '.': string += '.'
        return string

################################################################################

class Book(Reference):
    """
    A class for representing a complete book. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `authors`           A list of the authors of the reference
    `title`             The title of the reference
    `publisher`         The publisher of the book
    `address`           The address of the publisher (usually city and state/country)
    `volume`            The volume of the book
    `series`            The series the book belongs to
    `edition`           The edition of the book, as a string ordinal (e.g. ``'First'``)
    `year`              The year the reference was published (as a string)
    `url`               A DOI or other permalink to the reference
    =================== ========================================================

    """

    def __init__(self, authors=None, title='', publisher='', address='', volume='', series='', edition='', year='', url=''):
        Reference.__init__(self, authors=authors, title=title, year=year, url=url)
        self.publisher = publisher
        self.address = address
        self.volume = volume
        self.series = series
        self.edition = edition

    def __str__(self):
        """
        Return a string representation of the reference in reStructuredText
        format.
        """
        string = self.getAuthorString()
        if self.title != '':
            string += ' *%s.*' % (self.title)
        if self.edition != '':
            string += ' %s edition.' % (self.edition)
        if self.volume != '':
            string += ' Vol. %s.' % (self.volume)
        if self.address != '':
            string += ' %s:' % (self.address)
        if self.publisher != '':
            string += ' **%s**' % (self.publisher)
        if self.year != '':
            string += ' (%s)' % (self.year)
        return string + '.'

################################################################################

class Thesis(Reference):
    """
    A class for representing a graduate thesis. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `authors`           A list of the authors of the reference
    `title`             The title of the reference
    `degree`            ``'Ph.D.'`` or ``'Masters'``
    `school`            The name of the institution at which the thesis was written
    `year`              The year the reference was published (as a string)
    `url`               A DOI or other permalink to the reference
    =================== ========================================================

    """

    def __init__(self, authors=None, title='', degree='', school='', year='', url=''):
        Reference.__init__(self, authors=authors, title=title, year=year, url=url)
        self.degree = degree
        self.school = school
    
    def __str__(self):
        """
        Return a string representation of the reference in reStructuredText
        format.
        """
        string = self.getAuthorString()
        if self.title != '':
            string += ' "%s."' % (self.title)
        if self.degree != '':
            string += ' %s thesis.' % (self.degree)
        if self.school != '':
            string += ' ' % (self.school)
        if self.year != '':
            string += ' (%s)' % (self.year)
        if string[-1] != '.': string += '.'
        return string

################################################################################
