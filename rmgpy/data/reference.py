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
Contains classes and functions for working with the bibliographic information.
Currently there are three such classes:

* :class:`Article` - For articles in a journal, magazine, or other periodical

* :class:`Book` - For complete books

* :class:`Thesis` - For a graduate thesis

The above are all derived from the base :class:`Reference` class, which can
also be used if the reference does not fit into any of the above categories.
"""

import re


################################################################################

class Reference(object):
    """
    A base class for representing bibliographic information. The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `authors`           A list of the authors of the reference
    `title`             The title of the reference
    `year`              The year the reference was published (as a string)
    `doi`               A DOI link to the reference
    `url`               Any other link to the reference
    =================== ========================================================

    """

    def __init__(self, authors=None, title='', year='', doi='', url=''):
        self.authors = authors or []
        self.title = title
        self.year = year
        self.doi = doi
        self.url = url

    def __repr__(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = self.to_pretty_repr()
        string = re.sub(r'\(\n    ', '(', string)
        string = re.sub(r',\n    ', ', ', string)
        string = re.sub(r',\n\)', ')', string)
        string = re.sub(r' = ', '=', string)
        return string

    def __str__(self):
        """
        Return a string representation of the reference in reStructuredText
        format.
        """
        string = self.get_author_string()
        if self.title != '':
            string += f' *{self.title}*'
        if self.year != '':
            string += f' ({self.year})'
        if string and string[-1] != '.':
            string += '.'
        return string

    def to_pretty_repr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = 'Reference(\n'
        if len(self.authors) != 0:
            string += f"    authors = [{', '.join([f'{author!r}' for author in self.authors])}],\n"
        if self.title != '':
            string += f'    title = {self.title!r},\n'
        if self.year != '':
            string += f'    year = {self.year!r},\n'
        if self.doi != '':
            string += f'    doi = {self.doi!r},\n'
        if self.url != '':
            string += f'    url = {self.url!r},\n'
        return string + ')'

    def get_author_string(self):
        """
        Return a pretty, reStructuredText-formatted string of the authors.
        """
        authors = ''
        if self.authors is not None and len(self.authors) > 0:
            if len(self.authors) == 1:
                authors = f'{self.authors[0]}.'
            elif len(self.authors) == 2:
                authors = f'{self.authors[0]} and {self.authors[1]}.'
            elif self.authors[-1] == 'et al':
                authors = f"{', '.join(self.authors[:-1])} et al."
            else:
                authors = f"{', '.join(self.authors[:-1])}, and {self.authors[-1]}."
            # reStructuredText automatically interprets "A." et al as a 
            # numbered list; this suppresses that behavior
            if authors[1:3] == '. ':
                authors = authors[0:2] + '\ ' + authors[2:]
            # If the last author is of the form "Lastname, A. B.", this will
            # remove the extra period at the end of the sentence
            if authors[-2:] == '..':
                authors = authors[:-1]
        return authors


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
    `doi`               A DOI link to the reference
    `url`               Any other link to the reference
    =================== ========================================================

    """

    def __init__(self, authors=None, title='', journal='', volume='', number='', pages='', year='', doi='', url=''):
        Reference.__init__(self, authors=authors, title=title, year=year, doi=doi, url=url)
        self.journal = journal
        self.volume = volume
        self.number = number
        self.pages = pages

    def __str__(self):
        """
        Return a string representation of the reference in reStructuredText
        format.
        """
        string = self.get_author_string()
        if self.title != '':
            string += f' "{self.title}."'
        if self.journal != '':
            string += f' *{self.journal}*'
        if self.volume != '':
            string += f' **{self.volume}**'
        if self.number != '':
            string += f' ({self.number})'
        if self.pages != '':
            string += f', p. {self.pages}'
        if self.year != '':
            string += f' ({self.year})'
        if string and string[-1] != '.':
            string += '.'
        return string

    def to_pretty_repr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = 'Article(\n'
        if len(self.authors) != 0:
            string += f"    authors = [{', '.join([f'{author!r}' for author in self.authors])}],\n"
        if self.title != '':
            string += f'    title = {self.title!r},\n'
        if self.journal != '':
            string += f'    journal = {self.journal!r},\n'
        if self.volume != '':
            string += f'    volume = {self.volume!r},\n'
        if self.number != '':
            string += f'    number = {self.number!r},\n'
        if self.pages != '':
            string += f'    pages = {self.pages!r},\n'
        if self.year != '':
            string += f'    year = {self.year!r},\n'
        if self.doi != '':
            string += f'    doi = {self.doi!r},\n'
        if self.url != '':
            string += f'    url = {self.url!r},\n'
        return string + ')'


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
    `doi`               A DOI link to the reference
    `url`               Any other link to the reference
    =================== ========================================================

    """

    def __init__(self, authors=None, title='', publisher='', address='', volume='', series='', edition='', year='',
                 doi='', url=''):
        Reference.__init__(self, authors=authors, title=title, year=year, doi=doi, url=url)
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
        string = self.get_author_string()
        if self.title != '':
            string += f' *{self.title}.*'
        if self.edition != '':
            string += f' {self.edition} edition.'
        if self.volume != '':
            string += f' Vol. {self.volume}.'
        if self.address != '':
            string += f' {self.address}:'
        if self.publisher != '':
            string += f' **{self.publisher}**'
        if self.year != '':
            string += f' ({self.year})'
        return string + '.'

    def to_pretty_repr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = 'Book(\n'
        if len(self.authors) != 0:
            string += f"    authors = [{', '.join([f'{author!r}' for author in self.authors])}],\n"
        if self.title != '':
            string += f'    title = {self.title!r},\n'
        if self.publisher != '':
            string += f'    publisher = {self.publisher!r},\n'
        if self.address != '':
            string += f'    address = {self.address!r},\n'
        if self.volume != '':
            string += f'    volume = {self.volume!r},\n'
        if self.series != '':
            string += f'    series = {self.series!r},\n'
        if self.edition != '':
            string += f'    edition = {self.edition!r},\n'
        if self.year != '':
            string += f'    year = {self.year!r},\n'
        if self.doi != '':
            string += f'    doi = {self.doi!r},\n'
        if self.url != '':
            string += f'    url = {self.url!r},\n'
        return string + ')'


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
    `doi`               A DOI link to the reference
    `url`               Any other link to the reference
    =================== ========================================================

    """

    def __init__(self, authors=None, title='', degree='', school='', year='', doi='', url=''):
        Reference.__init__(self, authors=authors, title=title, year=year, doi=doi, url=url)
        self.degree = degree
        self.school = school

    def __str__(self):
        """
        Return a string representation of the reference in reStructuredText
        format.
        """
        string = self.get_author_string()
        if self.title != '':
            string += f' "{self.title}."'
        if self.degree != '':
            string += f' {self.degree} thesis.'
        if self.school != '':
            string += f' {self.school}'
        if self.year != '':
            string += f' ({self.year})'
        if string and string[-1] != '.':
            string += '.'
        return string

    def to_pretty_repr(self):
        """
        Return a string representation of the reference that can be used to
        reconstruct the object.
        """
        string = 'Thesis(\n'
        if len(self.authors) != 0:
            string += f"    authors = [{', '.join([f'{author!r}' for author in self.authors])}],\n"
        if self.title != '':
            string += f'    title = {self.title!r},\n'
        if self.degree != '':
            string += f'    degree = {self.degree!r},\n'
        if self.school != '':
            string += f'    school = {self.school!r},\n'
        if self.year != '':
            string += f'    year = {self.year!r},\n'
        if self.doi != '':
            string += f'    doi = {self.doi!r},\n'
        if self.url != '':
            string += f'    url = {self.url!r},\n'
        return string + ')'

################################################################################
