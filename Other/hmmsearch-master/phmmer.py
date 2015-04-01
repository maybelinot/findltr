#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib
import urllib2
from hmmerhit import HmmerHit
from fetchutil import readFasta, readAccession
import xml.etree.ElementTree as ET
import sys
import re
import glob
import os


class SmartRedirectHandler(urllib2.HTTPRedirectHandler):

    def http_error_301(
        self,
        req,
        fp,
        code,
        msg,
        headers,
        ):
        return headers

    def http_error_302(
        self,
        req,
        fp,
        code,
        msg,
        headers,
        ):
        return headers


class PhmmerSearch(object):

    """
....Search pdb sequence databases to extract matched PDB structure for a sequence using phmmer (http://hmmer.janelia.org/search/phmmer)

...."""

    def __init__(self, **kwargs):
        self.file = kwargs.get('file')
        self.name = kwargs.get('name')
        self.sequence = kwargs.get('sequence')
        self.source = kwargs.get('source')
        self.accession = kwargs.get('accession')
        self.cutoff = kwargs.get('evalue', 1e-5)
        self.length = kwargs.get('length')
        self.overwrite = kwargs.get('overwrite')
      
        self.species = kwargs.get('species')
        self.db = kwargs.get('db', 'pdb')
        self.tier = kwargs.get('tier')
        self.hits = []
        self.color = kwargs.get('color')

    def readFile(self, fileName):
        (self.name, self.sequence) = readFasta(fileName)
        if self.sequence:
            self.length = len(self.sequence)
            self.file = fileName
        (self.source, self.accession) = readAccession(self.file)
        self.phmmerSearch()
        return

    def removeRedundant(self):
        """
........Removal of overlapped PDB hits from hit list (self.hits)
........"""

        for i in range(len(self.hits)):
            for j in range(i + 1, len(self.hits)):
                one = self.hits[i]
                another = self.hits[j]

                # if 'one' has cover entire range of 'another' and its bitscore is higher than 'another', exclude 'another'
                # note that bitscore is stored as STRING so it should be typecasted as float before comparisions

                if one.start <= another.start and one.end \
                    >= another.end and float(one.bitscore) \
                    >= float(another.bitscore):
                    another.exclude = True

                # if 'another' has cover entire range of 'one' and its bitscore is higher than 'another', exclude 'one'

                if one.start >= another.start and one.end \
                    <= another.end and float(one.bitscore) \
                    <= float(another.bitscore):
                    one.exclude = True

                # if 'one' and 'another' is overlapped significant portions, (if overlapped length is 90% of length of domain)
                # exclude 'another'.....................

                if one.start <= another.start and one.end \
                    >= another.start or one.start >= another.start \
                    and one.start <= another.end:
                    overlaplength = min(one.end, another.end) \
                        - max(one.start, another.start)

                    # print one.start, one.end, another.start, another.end, overlaplength, one.end-one.start, one.end-one.start*0.8

                    if (one.end - one.start) * 0.9 < overlaplength:
                        another.exclude = True
                    overlaplength = 0

        return

    def phmmerParsing(self, result):

        stripchain = re.compile('(\S{4})_\S$')
        root = ET.fromstring(result)
        for child in root.iter('hits'):

            name = child.get('name')
            acc = child.get('acc')
            evalue = child.get('evalue')
            match = stripchain.match(name)
            if match:
                name = match.group(1).upper()

            for element in child.iter('domains'):
                if float(element.get('ievalue')) < float(self.cutoff):
                    cevalue = element.get('cevalue')
                    ievalue = element.get('ievalue')
                    start = element.get('alihmmfrom')
                    end = element.get('alihmmto')
                    bitscore = element.get('bitscore')
                    labellink = \
                        'http://www.rcsb.org/pdb/explore/explore.do?structureId={0}'.format(name)
                    hit = HmmerHit(
                        name=name,
                        acc=acc,
                        bitscore=bitscore,
                        labellink=labellink,
                        evalue=evalue,
                        ievalue=ievalue,
                        cevalue=cevalue,
                        start=start,
                        end=end,
                        )
                    hit.desc = 'feature'
                    hit.border = True
                    hit.startshow = False
                    hit.endshow = False
                    hit.label = True
                    hit.gradient = False
                    hit.tier = self.tier
                    self.hits.append(hit)
            self.removeRedundant()
            return

    def phmmerSearch(self):
        """
........search PDB hit using phmmer
........"""

        filename = self.file + '_pdb_phmmer.xml'
        if os.path.exists(filename) and not self.overwrite:
            print '{0} is already processed. Skipped.'.format(filename)
            f = open(filename)
            read = f.readlines()
            self.phmmerParsing(''.join(read))
            return True
        else:
            opener = urllib2.build_opener(SmartRedirectHandler())
            urllib2.install_opener(opener)
            print 'search {0} in remote database'.format(self.file)
            if not self.db in ['swissprot', 'pdb']:
                print 'invalid db. It should be either pdb or swissprot'
                sys.exit()

            # Strip chain name from PDB id

            parameters = {'seqdb': self.db, 'seq': self.sequence,
                          'algo': 'phmmer'}
            enc_params = urllib.urlencode(parameters)

            # post the seqrch request to the server

            request = \
                urllib2.Request('http://hmmer.janelia.org/search/phmmer'
                                , enc_params)

            # get the url where the results can be fetched from

            results_url = urllib2.urlopen(request).getheader('location')

            # Note that 'ali' parameter should be set as true otherwise service would not return alignment..

            res_params = {'output': 'xml', 'ali': 'true'}

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

                f = open(filename, 'w')
                f.write(result)
                f.close()
                self.phmmerParsing(result)
                return True
            else:
                print 'Failed to retrieve results'
                return False


if __name__ == '__main__':
    list = glob.glob('*.fasta')
    phmmerResults = []
    for item in list:
        print item
        phmmer = PhmmerSearch(db='pdb', tier=1, evalue=1e-5)
        phmmer.readFile(item)
        phmmerResults.append(phmmer)

        for pdb in phmmer.hits:
            if not pdb.exclude:
                print pdb.name, pdb.desc, pdb.start, pdb.end, \
                    pdb.bitscore, pdb.exclude

