#!/usr/bin/python
# -*- coding: utf-8 -*-

import re


def splitsentence(word):

    wordlist = re.split("(\W+)", word)

    wordcount = len(wordlist)
    splited = []
    print wordlist
    divpoint = wordcount // 2
    if wordcount > 2:
        splited.append(''.join(wordlist[:divpoint]))
        splited.append(''.join(wordlist[divpoint + 1:]))
    else:
        splited.append(word)
    return splited


if __name__ == '__main__':
    stupid = "This is really stupid, shit-isn't it?"
    print splitsentence(stupid)
