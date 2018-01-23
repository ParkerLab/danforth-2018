#!/usr/bin/env python


def load_pprint_dict(f):
    """Returns a dictionary, loaded from the representation in the passed file"""

    with open(f, 'r') as dsf:
        d = eval(dsf.read())

    return d
