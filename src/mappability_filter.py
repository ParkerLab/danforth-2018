#!/usr/bin/env python


def filter_bed(bed, whitelists=[], blacklists=[]):
    whitelist_template = 'intersectBed -a {} -b {} -f 1.0'
    whitelist_chain = ' | '.join([whitelist_template.format('stdin', x) for x in whitelists])

    blacklist_template = 'intersectBed -a {} -b {} -v'
    blacklist_chain = ' | '.join([blacklist_template.format('stdin', x) for x in blacklists])

    command = ' | '.join([x for x in ['cat {bed}'.format(**locals()), whitelist_chain, blacklist_chain] if x != ''])
    return command
