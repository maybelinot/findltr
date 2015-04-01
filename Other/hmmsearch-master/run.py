#!/usr/bin/python
# -*- coding: utf-8 -*-

from bottle import route, run, template


@route('/hello/:name')
def index(name):
    if int(name) < 6:

        list = u"{0}+{0}은 ".format(name)
        i = int(name)
        for k in range(i):
            list = list + u"귀요"
        list = list + u"미"
    else:
        list = u"그만해;;;;지겹다"
    return template('<b>{{name}}</b>!', name=list)

run(host='localhost', port=8080)
