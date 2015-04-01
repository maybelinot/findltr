#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib
import urllib2
import json
from hmmerhit import HmmerHit
from fetchutil import readFasta, readAccession
from svgdrawer import SVGDrawer
import re
import glob
import os


class SmartRedirectHandler(urllib2.HTTPRedirectHandler):

    def http_error_301(self, req, fp, code, msg, headers):
        return headers

    def http_error_302(self, req, fp, code, msg, headers):
        return headers


class miscAnnotation(object):

    def __init__(self, **kwargs):

        self.file = kwargs.get('file')
        self.name = kwargs.get('name')
        self.sequence = kwargs.get('sequence')
        self.source = kwargs.get('source')
        self.accession = kwargs.get('accession')
        self.length = kwargs.get('length')
        self.method = kwargs.get('method')
        self.tier = kwargs.get('tier')
        self.hits = []
        self.color = kwargs.get('color')
        self.overwrite = kwargs.get('overwrite')

    def readFile(self, fileName):
        (self.name, self.sequence) = readFasta(fileName)
        if self.sequence:
            self.length = len(self.sequence)
            self.file = fileName
        (self.source, self.accession) = readAccession(self.file)
        self.annotationReader()
        return

    def parseAnnotation(self, result):
        """
        Parse using ElementTree Modules (http://docs.python.org/2/library/xml.etree.elementtree.html)
        """
        extractRegex = re.compile('var sequence = (.*);')
        match = extractRegex.search(result)
        if match:
            jsonData = match.group(1)
            data = json.loads(jsonData)
            if 'motifs' in data:
                motifs = data['motifs']
                for motif in motifs:
                    desc = 'feature'
                    hit = HmmerHit(name=self.method, desc=desc,
                                   start=motif['start'], end=motif['end'])
                    hit.border = False
                    hit.startshow = False
                    hit.endshow = False
                    hit.label = False
                    hit.gradient = False
                    hit.color = self.color
                    hit.tier = self.tier
                    self.hits.append(hit)
        return

    def annotationReader(self):
        """
........Hmmer Misc Annotation Reader Class
........"""
        if self.sequence and len(self.sequence) > 0:
            filename = self.file + '_{0}.xml'.format(self.method)
            if os.path.exists(filename) and not self.overwrite:
                print '{0} is already processed. Skipped.'.format(filename)
                f = open(filename)
                read = f.readlines()
                self.parseAnnotation(''.join(read))
                return True
            else:
                opener = urllib2.build_opener(SmartRedirectHandler())
                urllib2.install_opener(opener)
                print 'Running {0} prediction of {1} using web service..'.format(self.method,
                        self.file)
                if not self.method in ['disorder', 'coils', 'phobius']:
                    print '{0} is not valid method. It should be one of these method : disorder, coils, phobius.'
                    print 'search will be carried out using disorder'
                    self.method = 'disorder'

                parameters = {'seq': self.sequence}
                enc_params = urllib.urlencode(parameters)

                # post the seqrch request to the server

                request = \
                    urllib2.Request('http://hmmer.janelia.org/annotation/{0}'.format(self.method),
                                    enc_params)

                # get the url where the results can be fetched from

                results_url = \
                    urllib2.urlopen(request).headers.dict['location']

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

                self.name = '{0} Prediction:'.format(self.method) + self.name
                if result:
                    f = open(filename, 'w')
                    f.write(result)
                    f.close()
                    self.parseAnnotation(result)
                return True
        else:
            print 'Failed to retrieve results'
            return False


if __name__ == '__main__':
    list = glob.glob('*.fasta')
    miscFeatures = {}
    miscResults = []
    for item in list:
        print item
        misc = miscAnnotation(method='coils', tier=1)
        misc.readFile(item)
        miscResults.append(misc)
    miscFeatures['coils']
    drawer = SVGDrawer(outputSVG=item + '_disorder.svg',
                       hmmerResults=miscResults)
    drawer.drawSVG()

