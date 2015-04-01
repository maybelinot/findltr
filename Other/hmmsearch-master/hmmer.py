#!/usr/bin/python
# -*- coding: utf-8 -*-

#
# Hmmer3 based Protein Domain Visualizer
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#

from __future__ import division
import os
import argparse
import glob
from hmmerhit import HmmerHit
from htmltable import HTMLTable
from inputfile import InputFile
from svgdrawer import SVGDrawer
from annotation import miscAnnotation
from uniprotannotation import UniprotAnnotation
from phmmer import PhmmerSearch
from psipred import PsipredAnnotation
from hmmscan import Hmmer
import SimpleHTTPServer
import SocketServer
import webbrowser

class SVGList(object):

    def __init__(self):

        self.header = '<br>'
        self.footer = '<br>'
        self.svgEmbeddingTemplate = \
            """
<div class="span12 scrollable domain" id="{0}">{1}</div><div class="span12 scrollable sequence">{2}</div>
"""
        self.svgEmbedContent = ''

    def svgContentFill(self, svgContent):
        fill = self.svgEmbeddingTemplate.format(*svgContent)
        self.svgEmbedContent = self.svgEmbedContent + fill

class HmmerScanRunner(object):
    '''
    Class for the Hmmer domain hits
    '''
    def __init__(self, **kwargs):
        self.inputFileName = kwargs.get('inputFileName')
        self.inputFile = None
        self.files = kwargs.get('files')
        self.db = kwargs.get('db')
        self.evalue = kwargs.get('evalue')
        self.local = kwargs.get('local')
        self.outputHTML = kwargs.get('outputHTML')
        self.outputSVG = kwargs.get('outputSVG')
        self.threshold = kwargs.get('threshold')
        self.scaleFactor = kwargs.get('scaleFactor', 1.5)
        self.hmmerResults = []
        self.titlemode = kwargs.get('titlemode', False)
        self.overwrite = kwargs.get('overwrite')

    def processHmmerResults(self):

        # inject information of input file into hmmerResults
        if self.inputFile:
            colors = [
                'aliceblue',
                'antiquewhite',
                'aqua',
                'aquamarine',
                'azure',
                'beige',
                'bisque',
                'black',
                'blanchedalmond',
                'blue',
                'blueviolet',
                'brown',
                'burlywood',
                'cadetblue',
                'chartreuse',
                'chocolate',
                'coral',
                'cornflowerblue',
                'cornsilk',
                'crimson',
                'cyan',
                'darkblue',
                'darkcyan',
                'darkgoldenrod',
                'darkgray',
                'darkgreen',
                'darkgrey',
                'darkkhaki',
                'darkmagenta',
                'darkolivegreen',
                'darkorange',
                'darkorchid',
                'darkred',
                'darksalmon',
                'darkseagreen',
                'darkslateblue',
                'darkslategray',
                'darkslategrey',
                'darkturquoise',
                'darkviolet',
                'deeppink',
                'deepskyblue',
                'dimgray',
                'dimgrey',
                'dodgerblue',
                'firebrick',
                'floralwhite',
                'forestgreen',
                'fuchsia',
                'gainsboro',
                'ghostwhite',
                'gold',
                'goldenrod',
                'gray',
                'green',
                'greenyellow',
                'grey',
                'honeydew',
                'hotpink',
                'indianred',
                'indigo',
                'ivory',
                'khaki',
                'lavender',
                'lavenderblush',
                'lawngreen',
                'lemonchiffon',
                'lightblue',
                'lightcoral',
                'lightcyan',
                'lightgoldenrodyellow',
                'lightgray',
                'lightgreen',
                'lightgrey',
                'lightpink',
                'lightsalmon',
                'lightseagreen',
                'lightskyblue',
                'lightslategray',
                'lightslategrey',
                'lightsteelblue',
                'lightyellow',
                'lime',
                'limegreen',
                'linen',
                'magenta',
                'maroon',
                'mediumaquamarine',
                'mediumblue',
                'mediumorchid',
                'mediumpurple',
                'mediumseagreen',
                'mediumslateblue',
                'mediumspringgreen',
                'mediumturquoise',
                'mediumvioletred',
                'midnightblue',
                'mintcream',
                'mistyrose',
                'moccasin',
                'navajowhite',
                'navy',
                'oldlace',
                'olive',
                'olivedrab',
                'orange',
                'orangered',
                'orchid',
                'palegoldenrod',
                'palegreen',
                'paleturquoise',
                'palevioletred',
                'papayawhip',
                'peachpuff',
                'peru',
                'pink',
                'plum',
                'powderblue',
                'purple',
                'red',
                'rosybrown',
                'royalblue',
                'saddlebrown',
                'salmon',
                'sandybrown',
                'seagreen',
                'seashell',
                'sienna',
                'silver',
                'skyblue',
                'slateblue',
                'slategray',
                'slategrey',
                'snow',
                'springgreen',
                'steelblue',
                'tan',
                'teal',
                'thistle',
                'tomato',
                'turquoise',
                'violet',
                'wheat',
                'white',
                'whitesmoke',
                'yellow',
                'yellowgreen',
                ]

            for i in range(len(self.inputFile.fileNames)):
                for hmmerResult in self.hmmerResults:
                    if hmmerResult.file == self.inputFile.fileNames[i]:
                        hmmerResult.name = \
                            self.inputFile.proteinNames[i]

            for domain in self.inputFile.domainDefinitions:
                for hmmerResult in self.hmmerResults:
                    if domain.proteinName == hmmerResult.name:
                        if domain.action == 'remove':
                            for (tier, hits) in \
                                hmmerResult.features.items():
                                for hit in hits:
                                    if hit.name == domain.domainName \
    and hit.start == domain.start and hit.end == domain.end:
                                        hit.exclude = True
                        if domain.action == 'new':
                            pseudoHit = \
                                HmmerHit(name=domain.domainName,
                                    start=domain.start, end=domain.end)
                            hmmerResult.features['domain'
                                    ].append(pseudoHit)

            for renameDef in self.inputFile.renameDefinitions:
                for hmmerResult in self.hmmerResults:
                    for (tier, hits) in hmmerResult.features.items():
                        for hit in hits:
                            if hit.name == renameDef.domainName:
                                hit.name = renameDef.newName

            for colorDef in self.inputFile.colorDefinitions:
                gradient = False
                border = True
                label = True
                number = True
                startshow = True
                endshow = True
                if colorDef.option:
                    options = colorDef.option.strip('+')
                    if 'gradient' in options:
                        gradient = True
                    if 'noborder' in options:
                        border = False
                    if 'nolabel' in options:
                        label = False
                    if 'nonumber' in options:
                        startshow = False
                        endshow = False
                    if 'nostart' in options:
                        startshow = False
                    if 'noend' in options:
                        endshow = False

                for hmmerResult in self.hmmerResults:
                    for (tier, hits) in hmmerResult.features.items():
                        for hit in hits:
                            if hit.name == colorDef.domainName \
                                and colorDef.color in colors:
                                hit.color = colorDef.color
                                hit.gradient = gradient
                                hit.border = border
                                hit.label = label
                                hit.number = number
                                hit.startshow = startshow
                                hit.endshow = endshow

        return

    def loadAnnotation(self):
        """
........Populate annotation data into hmmerResults.features dictionary.
........"""

        for hmmerResult in self.hmmerResults:
            disorder = miscAnnotation(method='disorder', tier=1,
                    color='grey',overwrite=self.overwrite)
            disorder.readFile(hmmerResult.file)
            hmmerResult.features['disorder'] = disorder.hits
            hmmerResult.tier[1] = 'Disorder'
            coils = miscAnnotation(method='coils', tier=2, color='green',
                                    overwrite=self.overwrite)
            coils.readFile(hmmerResult.file)
            hmmerResult.features['coils'] = coils.hits
            hmmerResult.tier[2] = 'Coiled-coil'
            uniprot = UniprotAnnotation(tier=3,overwrite=self.overwrite)
            uniprot.readFile(hmmerResult.file)
            hmmerResult.features['uniprot'] = uniprot.hits
            hmmerResult.tier[3] = 'UniProt'
            phmmer = PhmmerSearch(tier=4,overwrite=self.overwrite)
            phmmer.readFile(hmmerResult.file)
            hmmerResult.features['PDB'] = phmmer.hits
            hmmerResult.tier[4] = 'PDB'
            psipred = PsipredAnnotation(file=hmmerResult.file, tier=5)
            results = psipred.psipredReader()
            if results:
                hmmerResult.features['SS'] = psipred.hits
                hmmerResult.tier[5] = 'SS'

        return

    def run(self):

        if self.inputFileName:
            self.inputFile = InputFile(inputFileName=self.inputFileName)
            self.inputFile.readInputFile()
            files = self.inputFile.fileNames
            self.db = self.inputFile.db
            self.threshold = self.inputFile.threshold
            self.local = self.inputFile.local
            self.evalue = self.inputFile.evalue
            self.outputHTML = self.inputFile.outputHTML
            self.outputSVG = self.inputFile.outputSVG
        else:
            if not self.files:
                self.files = glob.glob('*.fasta')
            files = self.files

        if self.threshold:
            threshold = 'cut_ga'
        else:
            threshold = 'No'

        for file in files:
            hmmer = Hmmer(file=file, db=self.db, evalue=self.evalue,
                          threshold=threshold,overwrite=self.overwrite)
            if self.local:
                if hmmer.runLocal():
                    self.hmmerResults.append(hmmer)
            else:
                if hmmer.runRemote():
                    self.hmmerResults.append(hmmer)

        if len(self.hmmerResults) > 0:
            self.processHmmerResults()
            self.loadAnnotation()
            self.hmmerResults.sort(key=lambda x: x.name)
            draw = SVGDrawer(outputSVG=self.outputSVG,
                             scaleFactor=self.scaleFactor,
                             hmmerResults=self.hmmerResults,
                             outputHTML=self.outputHTML,
                             titlemode=self.titlemode)
            draw.drawSVG()
            (sequenceFileNames, sequenceContent) = draw.sequencesDraw()
            (svgFileNames, svgContent) = draw.drawMultiSVG()
            header = ['Accession', 'Name', 'Domain', 'length']
            table = HTMLTable(header=header)
            table.scriptcontent="""           
                    $('#listtable').dataTable({
                        "sDom": "<'row'<'span8'l>r>t<'row'<'span8'i><'span8'p>>",
                         "iDisplayLength": 50,
                         "aoColumnDefs": [
                         { "sWidth": "100px", "aTargets": [ 0 ] },
                         { "sWidth": "400px", "aTargets": [ 1 ] },
                         { "sWidth": "400px", "aTargets": [ 2 ] },
                         { "sWidth": "100px", "aTargets": [ 3 ] },
                                         
                         ]                                    
                        });
            """
            table.style ="""                  
                <style>
                table{
                    font-family: "Arial",Sans-Serif;
                    font-size: 12px;
                    margin: 40px;
                    width:1000px;
                    text-align: left;
                    border-collapse: collapse;  
                    }
                tr.conditionalRowColor
                {
                    border-left:2px;
                    border-left-color:red;
                    border-right:2px; 
                    border-right-color:red;
                    border-top:2px;
                    border-top-color:red;
                    border-bottom:2px 
                    border-botton-color:red;
                }
                    
                 td.conditionalRowColor
                {
                    background-color:#FFEEEE;
                }

                .scrollable {
                height: 100%;
                overflow: auto;
                }
                div.head {
                    width:800px;
                    font-family: Sans-Serif;
                    font-size: 14px;
                    border:3px solid #EEEEEE;
                    border-radius: 10px;
                    padding: 10px;
                    align :center;
                    background-color: #FFFFFF;
                    }
               div.dataTables_length label {
                    width: 460px;
                    float: left;
                    text-align: left;
                }
                 
                div.dataTables_length select {
                    width: 75px;
                }
                 
                div.dataTables_filter label {
                    float: right;
                    width: 460px;
                }
                 
                div.dataTables_info {
                    padding-top: 8px;
                }
                 
                div.dataTables_paginate {
                    float: right;
                    margin: 0;
                }
                 
                table {
                    clear: both;
                } 
                </style>
            """
            svgList = SVGList()

            for hmmerResult in self.hmmerResults:
                print hmmerResult.name, hmmerResult.accession, \
                    hmmerResult.source
                if hmmerResult.source == 'refseq':
                    linkAddress = \
                        "<a href='http://www.ncbi.nlm.nih.gov/protein/{0}'>{0}</a>"
                if hmmerResult.source == 'uniprot':
                    linkAddress = \
                        "<a href='http://www.uniprot.org/uniprot/{0}'>{0}</a>"
                if hmmerResult.source in ['refseq', 'uniprot']:
                    accessionLink = \
                        linkAddress.format(hmmerResult.accession)
                    nameLink = \
                        "<a href='#{1}'>{0}</a>".format(hmmerResult.name,
                            hmmerResult.accession)
                else:
                    accessionLink = ''
                    nameLink = ''

                domain = []
                #
                # Table formation
                # Todo : Template Engine
                #
                for hit in hmmerResult.features['domain']:
                    if hit.acc:
                        domainName = \
                            "<a href='http://pfam.sanger.ac.uk/family/{0}'>{1}</a>".format(hit.acc,
                                hit.name)
                    else:
                        domainName = hit.name
                    startEnd = '({0}-{1})'.format(hit.start, hit.end)
                    if hmmerResult.source == 'uniprot':
                        startEndLink = \
                            "<a href='http://www.uniprot.org/blast/?about={0}[{1}-{2}]'>{3}</a>".format(hmmerResult.accession,
                                hit.start, hit.end, startEnd)
                    else:
                        startEndLink = startEnd
                    domainName = domainName + startEndLink
                    domain.append(domainName)
                domains = ','.join(domain)
                length = hmmerResult.length
                table.tableContentFill([accessionLink, nameLink,
                        domains, length])
                svgList.svgContentFill([hmmerResult.accession,
                        svgContent[hmmerResult.name],
                        sequenceContent[hmmerResult.name]])
            svg = svgList.svgEmbedContent
            table.extra = '<br><div class="span12 scrollable">' + svg + '</div>'
            table.tableGenerate(self.outputHTML)
            PORT=8000
            linkFormat = 'http://localhost:'+str(PORT)+'/{0}'
            Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
            httpd = SocketServer.TCPServer(("", PORT), Handler)
            print "serving at port", PORT
            webbrowser.open(linkFormat.format(self.outputHTML))
            httpd.serve_forever()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--fasta',
        nargs='+',
        dest='files',
        default=[],
        help='Files to process',
        )
    parser.add_argument('-d', '--database', dest='db', default='pfam',
                        help='HMM database')
    parser.add_argument('-e', '--evalue_cutoff', dest='evalue',
                        default=1e-5, help='E-value cutoff')
    parser.add_argument(
        '-l',
        '--local',
        action='store_true',
        dest='local',
        default=False,
        help='run local Hmmer',
        )
    parser.add_argument('-o', '--outputhtml', dest='outputHTML',
                        default='output.html',
                        help='Output HTML filename')
    parser.add_argument('-s', '--outputsvg', dest='outputSVG',
                        default='output.svg', help='Output SVG filename'
                        )
    parser.add_argument(
        '-t',
        '--no_threshold',
        dest='threshold',
        action='store_false',
        default=True,
        help='Turn of Pfam gathering threshold. Enable to look up more weak(unreliable) domains'
        )
    parser.add_argument(
        '-w',
        '--overwrite',
        dest='overwrite',
        action='store_true',
        default=False,
        help='overwrite existing results'
        )
    parser.add_argument('-i', '--input_file', dest='inputFileName',
                        default='hmmer.INP',
                        help='Read configuration file')
    results = parser.parse_args()
    if not os.path.exists(results.inputFileName):
        results.inputFileName = None
    hmmerscan = HmmerScanRunner(
        files=results.files,
        inputFileName=results.inputFileName,
        db=results.db,
        evalue=results.evalue,
        local=results.local,
        outputHTML=results.outputHTML,
        outputSVG=results.outputSVG,
        threshold=results.threshold,
        overwrite=results.overwrite
        )
    hmmerscan.run()

