#!/usr/bin/python
# -*- coding: utf-8 -*-


class HmmerHit(object):

    """
....Hmmer domain hits data container

....name : Name of hit. It is used for the label display in svgdrawer.py
....acc : Accession #
....desc : Description
....evalue : evalue of hit
...."""

    def __init__(self, **kwargs):
        self.name = kwargs.get('name')
        self.acc = kwargs.get('acc')
        self.desc = kwargs.get('desc')
        self.evalue = kwargs.get('evalue')
        self.ievalue = kwargs.get('ievalue')
        self.cevalue = kwargs.get('cevalue')
        self.start = int(kwargs.get('start'))
        self.end = int(kwargs.get('end'))
        self.bitscore = kwargs.get('bitscore')
        self.labellink = kwargs.get('labellink')
        self.numberlink = kwargs.get('numberlink')
        self.tier = 0
        self.exclude = False
        self.color = None
        self.label = True
        self.border = True
        self.startshow = True
        self.endshow = True
        self.gradient = False


