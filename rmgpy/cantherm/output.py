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
This module contains helper functionality for writing CanTherm output files.
"""

import ast
    
################################################################################

class PrettifyVisitor(ast.NodeVisitor):
    """
    A class for traversing an abstract syntax tree to assemble a prettier
    version of the code used to create the tree. Used by the :func:`prettify`
    function.
    """
    
    def __init__(self, level=0, indent=4):
        self.string = ''
        self.level = level
        self.indent = indent
        
    def visit_Call(self, node):
        """
        Return a pretty representation of the class or function call 
        represented by `node`.
        """
        result = node.func.id + '(\n'
        
        keywords = []
        for keyword in node.keywords:
            keywords.append('{0}={1}'.format(keyword.arg, self.visit(keyword.value)))
        result = '{0}({1})'.format(node.func.id, ', '.join(keywords)) 
        
        if len(result) > 80:
            result = node.func.id + '(\n'
        
            self.level += 1
            for keyword in node.keywords:
                result += '{2}{0} = {1},\n'.format(keyword.arg, self.visit(keyword.value), ' ' * (self.level * self.indent))
            self.level -= 1

            result += ' ' * (self.level * self.indent) + ')'
            
        self.string = result
        
        return result
    
    def visit_List(self, node):
        """
        Return a pretty representation of the list represented by `node`.
        """
        if any([not isinstance(e, (ast.Str, ast.Num)) for e in node.elts]):
            # Split elements onto multiple lines
            result = '[\n'
            self.level += 1
            for element in node.elts:
                result += '{1}{0},\n'.format(self.visit(element), ' ' * (self.level * self.indent))
            self.level -= 1
            result += '{0}]'.format(' ' * (self.level * self.indent))
            return result
        else:
            # Keep elements on one line
            result = '[{0}]'.format(', '.join([self.visit(e) for e in node.elts]))
            self.string = result
            return result
    
    def visit_Tuple(self, node):
        """
        Return a pretty representation of the tuple represented by `node`.
        """
        # If the tuple represents a quantity, keep it on one line
        isQuantity = True
        if len(node.elts) == 0 or not isinstance(node.elts[0], (ast.Num,ast.List)) or (
            isinstance(node.elts[0], ast.List) and any([not isinstance(e, ast.Num) for e in node.elts[0].elts])):
            isQuantity = False
        elif len(node.elts) < 2 or not isinstance(node.elts[1], ast.Str):
            isQuantity = False
        
        if not isQuantity:
            # Split elements onto multiple lines
            result = '(\n'
            self.level += 1
            for element in node.elts:
                result += '{1}{0},\n'.format(self.visit(element), ' ' * (self.level * self.indent))
            self.level -= 1
            result += '{0})'.format(' ' * (self.level * self.indent))
            return result
        else:
            # Keep elements on one line
            result = '({0})'.format(', '.join([self.visit(e) for e in node.elts]))
            self.string = result
            return result
    
    def visit_Dict(self, node):
        """
        Return a pretty representation of the dict represented by `node`.
        """
        if any([not isinstance(e, (ast.Str, ast.Num)) for e in node.keys]) or any([not isinstance(e, (ast.Str, ast.Num)) for e in node.values]):
            # Split elements onto multiple lines
            result = '{\n'
            self.level += 1
            for key, value in zip(node.keys, node.values):
                result += '{2}{0}: {1},\n'.format(self.visit(key), self.visit(value), ' ' * (self.level * self.indent))
            self.level -= 1
            result += '{0}}}'.format(' ' * (self.level * self.indent))
            self.string = result
            return result
        else:
            # Keep elements on one line
            result = '{{{0}}}'.format(', '.join(['{0}: {1}'.format(self.visit(key), self.visit(value)) for key, value in zip(node.keys, node.values)]))
            self.string = result
            return result
    
    def visit_Str(self, node):
        """
        Return a pretty representation of the string represented by `node`.
        """
        result = repr(node.s)
        self.string = result
        return result
    
    def visit_Num(self, node):
        """
        Return a pretty representation of the number represented by `node`.
        """
        result = '{0:g}'.format(node.n)
        #result = repr(node.n)
        self.string = result
        return result

def prettify(string, indent=4):
    """
    Return a "pretty" version of the given `string`, representing a snippet of
    Python code such as a representation of an object or function. This 
    involves splitting of tuples, lists, and dicts (including parameter lists)
    onto multiple lines, indenting as appropriate for readability.
    """
    # Parse the node into an abstract syntax tree
    node = ast.parse(string)
    # Traverse the tree, assembling the pretty version of the string
    visitor = PrettifyVisitor(indent=indent)
    visitor.visit(node)
    # Return the pretty version of the string
    return visitor.string
