#!/usr/bin/python
# -*- coding: utf-8 -*-

#
# Simple HTML Table Generator
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#


class SVGList(object):

    def __init__(self, **kwargs):

        self.header = '<div>'
        self.footer = '</div>'
        self.svgEmbeddingTemplate = \
            """
<object id="{0}", data="{0}" type = "image/svg+xml"><br>
"""
        self.svgEmbedContent = ''

    def tableContentFill(self, tableContent):
        fill = self.tableContentTemplate.format(*tableContent)
        self.tableContent = self.tableContent + fill

    def tableGenerate(self, filename):
        table = self.htmlHeader + self.include + self.scriptheader \
            + self.scriptcontent + self.scriptfooter + self.style \
            + self.tableHeader + self.tableContent + self.tablefooter \
            + self.extra + self.htmlfooter
        f = open(filename, 'w')
        f.write(table)
        f.close()
        return ()


