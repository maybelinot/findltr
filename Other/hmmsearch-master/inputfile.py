#!/usr/bin/python
# -*- coding: utf-8 -*-
import re
import os
from fetchutil import fetchfromUniprot


class DomainDefinition(object):

    def __init__(self, **kwargs):
        self.proteinName = kwargs.get('proteinName')
        self.domainName = kwargs.get('domainName')
        self.accession = kwargs.get('acession')
        self.start = int(kwargs.get('start'))
        self.end = int(kwargs.get('end'))
        self.action = kwargs.get('action')


class ColorDefinition(object):

    def __init__(self, **kwargs):
        self.domainName = kwargs.get('domainName')
        self.accession = kwargs.get('acession')
        self.color = kwargs.get('color')
        self.option = kwargs.get('option')


class RenameDefinition(object):

    def __init__(self, **kwargs):
        self.domainName = kwargs.get('domainName')
        self.newName = kwargs.get('newName')


class InputFile(object):

    #
    # Input File (Hmmer.INP) processing class
    #

    def __init__(self, **kwargs):
        self.domainDefinitions = []
        self.colorDefinitions = []
        self.renameDefinitions = []
        self.inputFileName = kwargs.get('inputFileName')
        self.fileNames = []
        self.proteinNames = []
        self.db = 'pfam'
        self.threshold = False
        self.local = False
        self.evalue = '1e-5'
        self.outputHTML = 'output.html'
        self.outputSVG = 'output.svg'
        self.commands = {
            'FILE': 'file',
            'NEW': 'new',
            'RENAME': 'rename',
            'REMOVE': 'remove',
            'COLOR': 'color',
            'SIZE': 'size',
            'EVALUE': 'evalue',
            'DB': 'db',
            'SETUP': 'setup',
            }

    def setupLocal(self, subcontent):

        data = subcontent.split('\t')
        if len(data) > 0:
            self.local = True
            if os.path.exists(data[0]):
                self.db = data[0]
            else:
                print 'HMM database file {0} is not found. Remote Hmmscan will be perfomed.'.format(data[0])
                self.db = 'pfam'
                self.local = False
        else:
            print 'HMM database file {0} is not found. Remote Hmmscan will be perfomed.'.format(data[0])
            self.db = 'pfam'
            self.local = False
        return

    def setupRemote(self, subcontent):

        data = subcontent.split('\t')
        if len(data) > 0:
            self.local = False
            self.db = data[0]
        else:
            print 'HMM database name is not found. Remote Hmmscan will be perfomed using pfam', \
                format(data[0])
            self.db = 'pfam'
            self.local = False
        return

    def setupThreshold(self, subcontent):
        data = subcontent.split('\t')
        if len(data) > 0:
            if data[0].lower() == 'true':
                self.threshold = True
        else:
            self.threshold = False
        return

    def setupEvalue(self, subcontent):
        data = subcontent.split('\t')
        if len(data) > 0:
            self.evalue = data[0]
        else:
            self.evalue = '1e-5'
        return

    def setupOutputHTML(self, subcontent):
        data = subcontent.split('\t')
        if len(data) > 0:
            self.outputHTML = data[0]
        else:
            self.outputHTML = 'output.html'
        print data, len(data), data[0], self.outputHTML
        return

    def setupOutputSVG(self, subcontent):
        data = subcontent.split('\t')
        if len(data) > 0:
            self.outputSVG = data[0]
        else:
            self.outputSVG = 'output.svg'
        return

    def setup(self, content):
        subcommands = {
            'LOCAL': 'setupLocal',
            'REMOTE': 'setupRemote',
            'THRESHOLD': 'setupThreshold',
            'EVALUE': 'setupEvalue',
            'OUTPUTHTML': 'setupOutputHTML',
            'OUTPUTSVG': 'setupOutputSVG',
            }
        regex = re.compile("^(\S+)\s+(.*)$")
        match = regex.match(content)
        if match:
            keyword = match.group(1)
            subcontent = match.group(2)
            if keyword in subcommands:
                process = getattr(self, subcommands[keyword])
                process(subcontent)
        return

    def file(self, content):
        data = content.split('\t')
        uniprot = re.compile("uniprot:(\S+)")
        if len(data) > 1:
            name = data[0]
            fileName = data[1]
            match = uniprot.match(fileName)
            if match:
                fileName = fetchfromUniprot(match.group(1))
            if fileName:
                if os.path.exists(fileName):
                    self.fileNames.append(fileName)
                    self.proteinNames.append(name)
        return

    def new(self, content):
        data = content.split('\t')
        if len(data) > 1:
            proteinName = data[0]
            domainName = data[1]
            if data[2] < data[3]:
                start = data[2]
                end = data[3]
            else:
                start = data[3]
                end = data[2]

            if len(data) > 4:
                color = data[4]
                colorDef = ColorDefinition(domainName=domainName,
                        color=color)
                self.colorDefinitions.append(colorDef)
            if len(data) > 5:
                colorDef.option = data[5]

            domain = DomainDefinition(proteinName=proteinName,
                    domainName=domainName, start=start, end=end,
                    action='new')
            self.domainDefinitions.append(domain)
        return

    def rename(self, content):
        data = content.split('\t')
        if len(data) > 1:
            domainName = data[0]
            newName = data[1]

            if len(data) > 2:
                color = data[2]
                colorDef = ColorDefinition(domainName=newName,
                        color=color)
                self.colorDefinitions.append(colorDef)
            if len(data) > 3:
                colorDef.option = data[3]

            rename = RenameDefinition(domainName=domainName,
                    newName=newName)
            self.renameDefinitions.append(rename)
        return

    def remove(self, content):
        data = content.split('\t')
        if len(data) > 1:
            proteinName = data[0]
            domainName = data[1]
            if data[2] < data[3]:
                start = data[2]
                end = data[3]
            else:
                start = data[3]
                end = data[2]
            domain = DomainDefinition(proteinName=proteinName,
                    domainName=domainName, start=start, end=end,
                    action='remove')
            self.domainDefinitions.append(domain)
        return

    def color(self, content):
        data = content.split('\t')
        if len(data) > 1:
            domainName = data[0]
            color = data[1]
            colorDef = ColorDefinition(domainName=domainName,
                    color=color)
            self.colorDefinitions.append(colorDef)
            if len(data) > 2:
                colorDef.option = data[2]
        return

    def readInputFile(self):

        if os.path.exists(self.inputFileName):
            f = open(self.inputFileName, 'r')
            inputFileContent = f.readlines()
            regex = re.compile("^(\S+)\s+(.*)$")
            for line in inputFileContent:
                if line[0:1] != '#':
                    match = regex.match(line)
                    if match:
                        keyword = match.group(1)
                        content = match.group(2)
                        if keyword in self.commands:
                            process = getattr(self,
                                    self.commands[keyword])
                            process(content)
            return True
        else:
            return False


