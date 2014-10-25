#!/usr/bin/python
# -*- coding: utf-8 -*-

#
# Simple HTML Table Generator
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#


class HTMLTable(object):

    def __init__(self, **kwargs):

        self.htmlHeader = \
            """
<!DOCTYPE html>
<html>
	<head>
    	<title>List</title>
 		<meta charset='utf-8'>
	</head>
	<script type="text/javascript" charset="utf8" src="{0}"></script>
"""
        self.include = \
            """
		<!-- DataTables CSS -->
		<link rel="stylesheet" type="text/css" href="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/css/jquery.dataTables.css">
	    <!-- BootStrap CSS -->
        <link href="//netdna.bootstrapcdn.com/twitter-bootstrap/2.2.2/css/bootstrap-combined.min.css" rel="stylesheet">
		<!-- jQuery -->
		<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jQuery/jquery-1.8.2.min.js"></script>
	    
        <script src="//netdna.bootstrapcdn.com/twitter-bootstrap/2.2.2/js/bootstrap.min.js"></script>
		<!-- DataTables -->
		<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/jquery.dataTables.min.js"></script>
"""
        self.scriptheader = \
            """
		<script>
		$(document).ready(function(){		
"""
        self.scriptfooter = """ 		
		});
		</script>
"""
        self.scriptcontent = ''
        self.style = ""
        self.tablefooter = '</div></table>'
        self.htmlfooter = '</body></html>'
        self.tableContentTemplate = ''
        self.tableContentMarkTemplate =''
        self.tableContent = ''
        self.tableHeader = ''
        self.extra = ''
        self.header = kwargs.get('header')
        self.colNo = len(self.header)
        self.tableHeaderGenerate()

    def tableHeaderGenerate(self):
        self.tableHeader = "<body><div class=span12><table id='listtable' class='table table-striped table-bordered'><thead><tr>" \
            + ''.join(['<td>' + i + '</td>' for i in self.header]) \
            + '</tr></thead><tbody>'
        self.tableContentTemplate = '<tr>' + ''.join(['<td>{' + str(i)
                + '}</td>' for i in range(self.colNo)]) + '</tr>'
        self.tableContentMarkTemplate = "<tr class='conditionalRowColor'>"\
            + ''.join(["<td class='conditionalRowColor'>{" + str(i)\
            + '}</td>' for i in range(self.colNo)]) + '</tr>'

    def tableContentFill(self, tableContent):
        fill = self.tableContentTemplate.format(*tableContent)
        self.tableContent = self.tableContent + fill

    def tableContentMarkFill(self,tableContent):
        fill = self.tableContentMarkTemplate.format(*tableContent)     
        self.tableContent = self.tableContent + fill

    def tableGenerate(self, filename):
        table = self.htmlHeader + self.include + self.scriptheader \
            + self.scriptcontent + self.scriptfooter + self.style \
            + self.tableHeader + self.tableContent + self.tablefooter \
            + self.extra + self.htmlfooter
        f = open(filename, 'w')
        f.write(table)
        f.close()
        return ()


