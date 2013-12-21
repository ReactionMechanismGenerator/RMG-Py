#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 5, 2013

@author: rwest


'''

import hashlib
from os.path import normpath, walk, isdir, isfile, dirname, basename, \
    exists as path_exists, join as path_join

def path_checksum(paths):
    """
    Recursively calculates a checksum representing the contents of all files
    found with a sequence of file and/or directory paths.
    
    eg. path_checksum(['/tmp'])

    Based on post by David Moss  at
    http://code.activestate.com/recipes/576973-getting-the-sha-1-or-md5-hash-of-a-directory/
    """
    if not hasattr(paths, '__iter__'):
        raise TypeError('sequence or iterable expected not %r!' % type(paths))

    def _update_checksum(checksum, dirname, filenames):
        for filename in sorted(filenames):
            path = path_join(dirname, filename)
            if isfile(path):
                #print path
                fh = open(path, 'rb')
                while 1:
                    buf = fh.read(4096)
                    if not buf : break
                    checksum.update(buf)
                fh.close()

    chksum = hashlib.sha1()

    for path in sorted([normpath(f) for f in paths]):
        if path_exists(path):
            if isdir(path):
                walk(path, _update_checksum, chksum)
            elif isfile(path):
                _update_checksum(chksum, dirname(path), basename(path))

    return chksum.hexdigest()

if __name__ == '__main__':
    print path_checksum([r'/tmp', '/etc/hosts'])
