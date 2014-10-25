#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import glob
import re
from hmmerhit import HmmerHit
import xml.etree.cElementTree as ET


class PsipredAnnotation(object):

    """
....PSIPRED v3.0 Results Parser and Visualization class

...."""

    def __init__(self, **kwargs):

        self.file = kwargs.get('file')
        self.name = kwargs.get('name')
        self.sequence = kwargs.get('sequence')
        self.source = kwargs.get('source')
        self.accession = kwargs.get('accession')
        self.length = kwargs.get('length')
        self.tier = kwargs.get('tier')
        self.hits = []
        self.number = []
        self.residues = []
        self.prediction = []
        self.confidence = []

    def saveSVG(self, filename, doc):
        """
........Save doc (ElementTreeDoc) content as a svg file.
........"""

        f = open(filename, 'w')
        f.write('<?xml version="1.0" standalone="no"?>\n')
        f.write('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n')
        f.write('"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')
        f.write('\n'.join(ET.tostringlist(doc)))
        f.close()
        return ()

    def secondaryDraw(self):
        """
........Draw Secondary Structure
........"""

        y = 70
        leftMargin = 200
        fontSize = 16
        columnWidth = 60
        canvasWidth = 1200
        yDelta = 30
        predictText = {'C': ' ', 'H': u"█", 'E': u"➜"}
        filenamebaseRegex = re.compile("(\S+).fasta")
        match = filenamebaseRegex.match(self.file)
        if match:
            svgfilename = match.group(1) + '.svg'
        else:
            svgfilename = 'output.svg'
        print self.file
        canvasHeight = int((self.length / float(columnWidth) + 2)
                           * yDelta)

        # doc is elementTree container for svg

        doc = ET.Element('svg', width=str(canvasWidth),
                         height=str(canvasHeight), version='1.2',
                         xmlns='http://www.w3.org/2000/svg')
        doc.attrib['xmlns:xlink'] = 'http://www.w3.org/1999/xlink'

        cursor = 0
        while cursor < self.length:
            numbers = ET.Element('text', x=str(leftMargin - fontSize),
                                 y=str(y), fill='black',
                                 style='font-family:Courier;font-size:'
                                 + str(fontSize)
                                 + 'px;text-anchor:end;dominant-baseline:middle'
                                 )
            numbers.text = str(self.number[cursor])
            doc.append(numbers)
            text = self.sequence[cursor:cursor + columnWidth]
            predict = self.prediction[cursor:cursor + columnWidth]
            aminoacids = ET.Element('text', x=str(leftMargin),
                                    y=str(y), fill='black',
                                    style='font-family:Courier;font-size:'
                                     + str(fontSize)
                                    + 'px;text-anchor:left;dominant-baseline:middle'
                                    )
            aminoacids.text = text
            doc.append(aminoacids)

            for i in range(len(predict)):
                if predict[i] == 'C':
                    color = 'black'
                if predict[i] == 'H':
                    color = 'green'
                if predict[i] == 'E':
                    color = 'red'

                prediction = ET.Element('text', x=str(leftMargin + i
                        * fontSize * 0.625), y=str(y + fontSize),
                        fill=color,
                        style='font-family:Courier;font-size:'
                        + str(fontSize)
                        + 'px;text-anchor:left;dominant-baseline:middle'
                        )
                prediction.text = predictText[predict[i]]
                doc.append(prediction)

            if cursor + columnWidth < self.length:
                endnumbers = ET.Element('text', x=str(leftMargin
                        + columnWidth * fontSize * 0.68), y=str(y),
                        fill='black',
                        style='font-family:Courier;font-size:'
                        + str(fontSize)
                        + 'px;text-anchor:end;dominant-baseline:middle')
                endnumbers.text = str(int(self.number[cursor])
                        + columnWidth - 1)
                doc.append(endnumbers)
            cursor += columnWidth
            y += yDelta

        predict = self.prediction[cursor:self.length]
        text = self.sequence[cursor:self.length]
        aminoacids = ET.Element('text', x=str(leftMargin), y=str(y),
                                fill='black',
                                style='font-family:Courier;font-size:'
                                + str(fontSize)
                                + 'px;text-anchor:left;dominant-baseline:middle'
                                )
        aminoacids.text = text
        prediction = ET.Element('text', x=str(leftMargin), y=str(y
                                + fontSize), fill='black',
                                style='font-family:Courier;font-size:'
                                + str(fontSize)
                                + 'px;text-anchor:left;dominant-baseline:middle'
                                )
        prediction.text = ''.join(predict)
        doc.append(aminoacids)
        doc.append(prediction)
        self.saveSVG(svgfilename, doc)

        return

    def psipredReader(self):
        """
........PsiPred output reader
........"""

        state = {'C': 'Coil', 'H': 'Helix', 'E': 'Sheet'}
        filenamebaseRegex = re.compile("(\S+).fasta")
        match = filenamebaseRegex.match(self.file)
        if match:
            filename = match.group(1) + '.ss'
        else:
            return False

        if os.path.exists(filename):
            print 'Reading PSIPRED File {0}...'.format(filename)
            f = open(filename, 'r')
            data = f.readlines()
            f.close()
            for line in data:
                content = line.split()
                self.number.append(content[0])
                self.residues.append(content[1])
                self.prediction.append(content[2])
                if content[2] == 'C':
                    self.confidence.append(content[3])
                if content[2] == 'H':
                    self.confidence.append(content[4])
                if content[2] == 'E':
                    self.confidence.append(content[5])

            currentstate = self.prediction[0]
            start = self.number[0]
            end = self.number[0]
            for i in range(1, len(self.number)):
                if currentstate != self.prediction[i]:
                    self.addhit(start, state, end, currentstate)

                    start = self.number[i]
                currentstate = self.prediction[i]
                end = self.number[i]

            self.addhit(start, state, end, currentstate)
            self.sequence = ''.join(self.residues)
            self.length = len(self.sequence)
            return True
        else:
            return False

    def addhit(
        self,
        start,
        state,
        end,
        currentstate,
        ):
        if currentstate != 'C':
            hit = HmmerHit(start=start, end=end, desc='feature',
                           name=state[currentstate])
            hit.tier = self.tier
            hit.label = False
            hit.border = False
            hit.startshow = False
            hit.endshow = False
            hit.gradient = False
            if currentstate == 'H':
                color = 'green'
            else:
                color = 'red'
            hit.color = color
            self.hits.append(hit)


if __name__ == '__main__':

    list = glob.glob('*.fasta')
    uniprotResults = []
    for item in list:
        psipred = PsipredAnnotation(file=item)
        success = psipred.psipredReader()
        if success:
            for hit in psipred.hits:
                print hit.start, hit.end, hit.name
            psipred.secondaryDraw()

