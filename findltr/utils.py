#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: maybelinot
# @Email: edik.trott@yandex.ru
# @Date:   2015-09-12 16:06:18
# @Last Modified by:   maybelinot
# @Last Modified time: 2015-09-12 20:14:58

from __future__ import unicode_literals, absolute_import

import logging
import os
import subprocess
import sys
import time

# EXTERNALLY INSTALLED
from BCBio import GFF
from Bio import SeqIO, Seq, SeqRecord, SeqFeature
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
import yaml

# Load logging before anything else
logging.basicConfig(format='>> %(message)s')
logr = logging.getLogger('findltr')

''' Load the config file so modules can import and reuse '''
CONFIG_FILE = os.path.expanduser('~/.findltr')
if os.path.exists(CONFIG_FILE):
    with open(CONFIG_FILE) as _:
        config = yaml.load(_)
else:
    config = {}


def export_gff(seq, young_lcp):

    unique_name = 'rec_%s.gff' % time.time()
    logr.info('Found LTRs are saved in ' + unique_name)

    records = []
    # fix name to chrN
    gff = SeqRecord.SeqRecord(Seq.Seq(seq), "seq0")
    top_feature = []
    for idx, item in enumerate(young_lcp):
        seq1 = SeqRecord.SeqRecord(
            Seq.Seq(seq[item[0][0]:item[0][1]]), id="seq1")
        seq2 = SeqRecord.SeqRecord(
            Seq.Seq(seq[item[1][0]:item[1][1]]), id="seq2")

        with open("/tmp/seq1.fasta", "w") as query:
            SeqIO.write(seq1, query, "fasta")
        with open("/tmp/seq2.fasta", "w") as subject:
            SeqIO.write(seq2, subject, "fasta")

        blast_output = NcbiblastnCommandline(
            query="/tmp/seq1.fasta", subject="/tmp/seq2.fasta", outfmt=5)()[0]
        blast_result_record = NCBIXML.read(StringIO(unicode(blast_output, 'utf-8')))
        identity = 0
        for alignment in blast_result_record.alignments:
            for hsp in alignment.hsps:
                identity = max(
                    hsp.identities / float(hsp.align_length) * 100.0, identity)
        identity = "%0.2f" % identity
        # cut zeros tail
        # identity = identity.rstrip("0")
        # identity = identity.rstrip(".")

        sub_qualifiers = {"source": "ltrfind", "ID": "UnknownLTR_" +
                          str(idx + 1), "ltr_similarity": identity}
        top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[0][0], item[1][1]),
                                                 type="Retrotransposon", strand=0, qualifiers=sub_qualifiers))

        sub_qualifiers_ltrs = {"source": "ltrfind", "Parent": "UnknownLTR_" +
                          str(idx + 1)}

        top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[0][0], item[0][1]),
                                                 type="long_terminal_repeat", strand=0, qualifiers=sub_qualifiers_ltrs))

        top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[1][0], item[1][1]),
                                                 type="long_terminal_repeat", strand=0, qualifiers=sub_qualifiers_ltrs))
    gff.features = top_feature
    with open(unique_name, "w") as out_handle:
        GFF.write([gff], out_handle)


class StandardArgs(object):

    '''
    FIXME: DOCS...
    '''
    uid = None
    password = None
    base_url = None

    def __init__(self, args=None, config=None):
        '''
        FIXME: DOCS...
        '''
        args = args or {}
        config = config or {}
        assert isinstance(args, dict)
        assert isinstance(config, dict)
        self._args, self._config = args, config

    def get(self, key, default=None, **kwargs):
        user_value = self._args.get(key) or self._config.get(key)
        if hasattr(default, '__call__'):
            value = user_value or default(**kwargs)
        else:
            value = user_value or default
        return value


def run(cmd):
    cmd = cmd if isinstance(cmd, list) else cmd.split()
    try:
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as error:
        logr.error("'{0}' failed: {1}".format(cmd, error))
        raise
    output, errors = process.communicate()
    if process.returncode != 0 or errors:
        if output:
            logr.error(output)
        if errors:
            logr.error(errors)
        sys.exit(process.returncode)
    return output, errors
