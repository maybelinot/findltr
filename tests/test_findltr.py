#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: maybelinot
# @Email: edik.trott@yandex.ru
# @Date:   2015-09-16 21:15:46
# @Last Modified by:   maybelinot
# @Last Modified time: 2015-09-16 21:18:33

from __future__ import unicode_literals
import pytest


def test_global_import():
    from findltr import algorithms
    from findltr import grouping
    from findltr import utils
