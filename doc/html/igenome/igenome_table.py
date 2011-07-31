#!/usr/bin/env python

import sys, os
import string
from ftplib import FTP

use_message = '''
'''

ftp_address = 'ftp.illumina.com'
user_id = 'igenome'
password = 'G3nom3s4u'

visible_organisms = \
    ['Homo_sapiens',
     'Mus_musculus',
     'Rattus_norvegicus',
     'Bos_taurus',
     'Canis_familiaris',
     'Gallus_gallus',
     'Drosophila_melanogaster',
     'Arabidopsis_thaliana',
     'Caenorhabditis_elegans',
     'Saccharomyces_cerevisiae']

def generate_table():
    ftp = FTP(ftp_address)

    print >> sys.stderr, "Logging into %s" % ftp_address
    print >> sys.stderr, ftp.login(user_id, password)

    organisms = {}
    organisms_index = {}

    def get_list(size_and_date = False):
        result = []
        lines = []
        def get_string(str):
            if len(str) > 0:
                lines.append(str)
                
        ftp.retrlines('LIST', get_string)
        
        for line in lines:
            elements = line.split()
            flag, organism, size, date = elements[0], elements[-1], elements[4], ' '.join(elements[-4:-1])
            result.append([flag, organism])

            if size_and_date:
                result[-1] += [size, date]

        return result

    files = get_list()
    for file in files:
        flag, name = file
        if flag[0] == 'd':
            organisms[name] = []
            organisms_index[name] = []

    for name in visible_organisms:
        print >> sys.stderr, name

        count = 0
        ftp.cwd("/%s" % name)
        sources = get_list()
        for flag, source in sources:
            if flag[0] != 'd':
                pass
            ftp.cwd("/%s/%s" % (name, source))
            versions = get_list()
            for flag, version in versions:
                if flag[0] != 'd':
                    next

                ftp.cwd("/%s/%s/%s" % (name, source, version))
                size_and_date = True
                files = get_list(True)
                file_name = ""
                for flag, file, size, date in files:
                    if string.find(file, 'tar.gz') != -1:
                        file_name = file
                        break

                size = int(size) / (1024 * 1024)
                size = "%d MB" % size
                    
                address = r"ftp://%s:%s@%s/%s/%s/%s/%s" % \
                    (user_id, password, ftp_address, name, source, version, file_name)
                organisms[name].append([source, version, address, size, date])
                count += 1
            organisms_index[name].append(count)
            
    ftp.close()

    html_file = open('igenome_table.html', 'w')
    print >> html_file, r'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">'
    print >> html_file, r'<html xmlns="http://www.w3.org/1999/xhtml">'
    print >> html_file, r'<head>'
    print >> html_file, r'<title>TopHat :: Center for Bioinformatics and Computational Biology</title>'
    print >> html_file, r'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />'
    print >> html_file, r'<link rel="stylesheet" type="text/css" href="css/style.css" media="screen" />'
    print >> html_file, r'</head>'
    print >> html_file, r'<body>'
    print >> html_file, r'<div id="leftside">'
    print >> html_file, r'<TABLE CELLPADDING=3 BORDER="1">'
    print >> html_file, r'<TR>'
    print >> html_file, r'<TH ALIGN="LEFT">Organism</TH>'
    print >> html_file, r'<TH ALIGN="LEFT">Data source</TH>'
    print >> html_file, r'<TH ALIGN="LEFT">Version</TH>'
    print >> html_file, r'<TH ALIGN="LEFT">Size</TH>'
    print >> html_file, r'<TH ALIGN="LEFT">Last Modified</TH>'
    print >> html_file, r'</TR>'

    for name in visible_organisms:
        organism = organisms[name]
        organism_index = organisms_index[name]

        prev_i = 0
        for i in organism_index:
            for j in range(prev_i, i):
                version = organism[j]
                print >> html_file, r'<TR>'

                if j == 0:
                    print >> html_file, r'<TD ALIGN="LEFT" rowspan=%d>%s</TD>' % \
                        (organism_index[-1], string.replace(name, '_', ' '))

                if j == prev_i:
                    print >> html_file, r'<TD ALIGN="LEFT" rowspan=%d>%s</TD>' % (i - prev_i, version[0])

                print >> html_file, r'<TD ALIGN="LEFT"><A HREF="%s">%s</A></TD>' % (version[2], version[1])
                print >> html_file, r'<TD ALIGN="RIGHT">%s</TD>' % version[3]
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD>' % version[4]
                print >> html_file, r'</TR>'
            
            prev_i = i
    
    print >> html_file, r'</TABLE>'
    print >> html_file, r'</div>'
    print >> html_file, r'</body>'
    print >> html_file, r'</html>'

    html_file.close()

if __name__ == "__main__":
    sys.exit(generate_table())
