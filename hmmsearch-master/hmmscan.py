#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Hmmerscan wrapper
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#
from __future__ import division
import urllib
import urllib2
from urllib2 import HTTPError
import os
import xml.etree.cElementTree as ET
import sys
import subprocess
from hmmerhit import HmmerHit
from abstractsequenceobject import AbstractSequenceObject
# install a custom handler to prevent following of redirects automatically.

class SmartRedirectHandler(urllib2.HTTPRedirectHandler):

    def http_error_302(self, req, fp, code, msg, headers):
        return headers

class Hmmer(AbstractSequenceObject):

    '''
    Hmmscan wrapper class
    '''

    def __init__(self, **kwargs):

        super(Hmmer, self).__init__(**kwargs)
        self.cutoff = kwargs.get('evalue',0.01)
        self.localHmmDB = kwargs.get('localHmmDB')
        self.threshold = kwargs.get('threshold')
        self.overwrite = kwargs.get('overwrite')
        self.excludeRedundancy = kwargs.get('excludeRedundancy',True)
        hits = []
        self.features['domain'] = hits
        self.tier = {}
        self.tier[0] = 'Pfam'

    def exclude(self):

        #
        # Exclude overlapped domain. If two domains are overlapped, the one have higher bitscore will be retained
        #
        hits = self.features['domain']

        if len(hits) > 1:
            for i in range(len(hits) - 1):
                for j in range(i + 1, len(hits)):
                    if not hits[i].exclude and not hits[j].exclude:
                        if hits[i].end > hits[j].start and hits[i].end < hits[j].end:
                            if hits[i].score > hits[j].score:
                                hits[j].exclude = True
                            else:
                                hits[i].exclude =True                  

                        if hits[j].end > hits[i].start and hits[j].end < hits[i].end:
                            if hits[i].score > hits[j].score :
                                hits[j].exclude = True
                            else:
                                hits[i].exclude = True

    def runRemote(self):

        '''
        Using Hmmscan in Hmmer3 web service, find locations of domains in the Fasta sequence and store into class
        '''

        if os.path.exists(self.file + '.xml') and not self.overwrite:
            print '{0} is already processed. Skipped.'.format(self.file
                    + '.xml')
            f = open(self.file + '.xml')
            read = f.readlines()
            self.parseHmmerScanXML(''.join(read))
            return True
        else:
            opener = urllib2.build_opener(SmartRedirectHandler())
            urllib2.install_opener(opener)
            print 'Running Hmmscan of {0} using web service..'.format(self.file)
            if not self.db in ['pfam', 'superfamily', 'tigrfam',
                               'gene3d']:
                print '{0} is not valid db. It should be pfam, superfamily, tigrfam or gene3d.'
                print 'search will be carried out with pfam'
                self.db = 'pfam'

            parameters = {'hmmdb': self.db,
                          'threshold': self.threshold,
                          'seq': self.sequence}
            enc_params = urllib.urlencode(parameters)

            try:

                # post the seqrch request to the server

                request = \
                    urllib2.Request('http://hmmer.janelia.org/search/hmmscan'
                                    , enc_params)
                # get the url where the results can be fetched from
                results_url = \
                    urllib2.urlopen(request).getheader('location')
                # modify the range, format and presence of alignments in your results here
                res_params = {'output': 'xml'}

                # add the parameters to your request for the results

                enc_res_params = urllib.urlencode(res_params)
                modified_res_url = results_url + '?' + enc_res_params

                # send a GET request to the server

                results_request = urllib2.Request(modified_res_url)
                data = urllib2.urlopen(results_request)

                # print out the results

                result = data.read()
                if result:

                    #
                    # Parse using ElementTree Modules (http://docs.python.org/2/library/xml.etree.elementtree.html)
                    #

                    f = open(self.file + '.xml', 'w')
                    f.write(result)
                    f.close()

                    self.parseHmmerScanXML(result)
                    return True
                else:
                    print 'Failed to retrieve results'
                    return False
            except HTTPError:
                print 'Hmmscan error'
                return False

    def parseHmmerScanXML(self, result):
        '''
        Parse Hmmerscan XML output into 
        '''

        root = ET.fromstring(result)
        for child in root.iter('hits'):
            name = child.get('name')
            desc = child.get('desc')
            acc = child.get('acc')
            evalue = child.get('evalue')

            for element in child.iter('domains'):

                if float(element.get('ievalue')) < float(self.cutoff):
                    cevalue = element.get('cevalue')
                    ievalue = element.get('ievalue')
                    start = element.get('iali')
                    end = element.get('jali')
                    bitscore = element.get('bitscore')
                    labellink = \
                        'http://pfam.sanger.ac.uk/family/{0}'.format(acc)

                    if self.source == 'uniprot':
                        numberlink = \
                            'http://www.uniprot.org/blast/?about={0}[{1}-{2}]'.format(self.accession,
                                start, end)
                    else:
                        numberlink = ""
                    hit = HmmerHit(
                        name=name,
                        desc=desc,
                        acc=acc,
                        bitscore=bitscore,
                        labellink=labellink,
                        numberlink=numberlink,
                        evalue=evalue,
                        ievalue=ievalue,
                        cevalue=cevalue,
                        start=start,
                        end=end,
                        )
                    self.features['domain'].append(hit)

        if self.excludeRedundancy:
            self.exclude()
        return ()

    def runLocal(self):

        '''
        Using Hmmscan in locally installed Hmmer3 packages, find locations of domains in the Fasta sequence and store
        '''

        if os.path.exists(self.db):
            outputFile = self.file + '.tb'
            hmmeroutFile = self.file + '.hmmer'
            if not os.path.exists(outputFile) or self.overwrite:
                print 'Running Hmmscan of {0} using local hmmer3'.format(self.file)
                #
                # hmmscan --domtbout outputFile --cut_ga, DBfile, queryFile
                #
                if self.threshold == 'cut_ga':
                    p = subprocess.Popen([
                        'hmmscan',
                        '--domtblout',
                        outputFile,
                        '--cut_ga',
                        self.db,
                        self.file,
                        ], stdout=subprocess.PIPE)
                else:
                    p = subprocess.Popen(['hmmscan', '--domtblout',
                            outputFile, self.db, self.file],
                            stdout=subprocess.PIPE)

                p_stdout = p.stdout.read()
                fw = open(hmmeroutFile, 'w')
                fw.write(p_stdout)
                fw.close()

            f = open(outputFile, 'r')
            data = f.readlines()
            f.close()

            #
            # Parse table formatted domain datafile generated by '--domtblout'.
            #
            descriptionLocation = 192
            for line in data:
                if line[0:1] != '#':
                    splited = line.split()
                    if float(splited[12]) < float(self.cutoff):
                        cevalue = splited[11]
                        ievalue = splited[12]
                        start = splited[17]
                        end = splited[18]
                        bitscore = splited[7]
                        acc = splited[1]
                        desc = line[descriptionLocation:]
                        evalue = splited[6]
                        name = splited[0]
                        query=splited[3]
                        labellink = \
                            'http://pfam.sanger.ac.uk/family/{0}'.format(acc)
                        if self.source == 'uniprot':
                            positionlink = \
                                'http://www.uniprot.org/blast/?about={0}[{1}-{2}]'.format(self.accession,
                                    start, end)
                            hit = HmmerHit(
                                name=name,
                                desc=desc,
                                acc=acc,
                                bitscore=bitscore,
                                labellink=labellink,
                                numberlink=positionlink,
                                evalue=evalue,
                                ievalue=ievalue,
                                cevalue=cevalue,
                                start=start,
                                end=end,
                                query=query
                                )
                        else:
                            hit = HmmerHit(
                                name=name,
                                desc=desc,
                                acc=acc,
                                bitscore=bitscore,
                                labellink=labellink,
                                evalue=evalue,
                                ievalue=ievalue,
                                cevalue=cevalue,
                                start=start,
                                end=end,
                                )
                        self.features['domain'].append(hit)
                else:
                    if line[0:8]=="# target":
                        descriptionLocation = line.find('description of target')
            if self.excludeRedundancy:
                self.exclude()
            return True
        else:
            print '{0} file is not exist'.format(self.db)
            sys.exit()
            return False


