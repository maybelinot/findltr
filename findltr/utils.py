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


def export_gff(seq, young_lcp, outputfile):

    gff_output = outputfile or 'rec_%s.gff3' % time.time()
    logr.info('Found LTRs are saved in ' + gff_output)

    records = []
    # fix name to chrN based on input seq
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
        blast_result_record = NCBIXML.read(StringIO(unicode(blast_output, "utf-8")))
        identity = 0
        for alignment in blast_result_record.alignments:
            for hsp in alignment.hsps:
                identity = max(
                    hsp.identities / float(hsp.align_length) * 100.0, identity)
        identity = "%0.2f" % identity
        # cut zeros tail
        # identity = identity.rstrip("0")
        # identity = identity.rstrip(".")
        # improve seqfeatures appending

        sub_qualifiers_region = {"source": "ltrfind",
                            "ID": "repeat_region" + str(idx + 1)}
        top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[0][0] - 4, item[1][1] + 4),
                                 type="repeat_region", strand=0, qualifiers=sub_qualifiers_region))

        sub_qualifiers_target_site = {"source": "ltrfind",
                            "Parent": "repeat_region" + str(idx + 1)}
        top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[0][0] - 4, item[0][0]),
                                 type="target_site_duplication", strand=0, qualifiers=sub_qualifiers_target_site))
        sub_qualifiers = {"source": "ltrfind",
                            "ID": "LTR_retrotransposon" + str(idx + 1),
                            "Parent": "repeat_region" + str(idx + 1),
                            "ltr_similarity": identity,
                            "seq_number": "0"}
        top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[0][0], item[1][1]),
                                 type="LTR_retrotransposon", strand=0, qualifiers=sub_qualifiers))

        sub_qualifiers_ltrs = {"source": "ltrfind", "Parent": "LTR_retrotransposon" +
                          str(idx + 1)}

        top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[0][0], item[0][1]),
                                 type="long_terminal_repeat", strand=0, qualifiers=sub_qualifiers_ltrs))

        top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[1][0], item[1][1]),
                                 type="long_terminal_repeat", strand=0, qualifiers=sub_qualifiers_ltrs))

        top_feature.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(item[1][1], item[1][1] + 4),
                                 type="target_site_duplication", strand=0, qualifiers=sub_qualifiers_target_site))
    gff.features = top_feature
    # track name='findltr' description='findltr Supplied Track'

    with open(gff_output, "w") as out_handle:
        GFF.write([gff], out_handle)


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
