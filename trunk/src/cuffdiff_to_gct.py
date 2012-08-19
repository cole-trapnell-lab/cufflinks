#!/usr/bin/env python
# encoding: utf-8
"""
cuffdiff_to_gct.py

Created by Cole Trapnell on 2010-08-30.
Copyright (c) 2010 Cole Trapnell. All rights reserved.
"""

import sys
import getopt
import math

help_message = '''
This script transforms a Cuffdiff FPKM tracking file to an IGV-compatible GCT expression file

Usage:
    cuffdiff_to_gct.py <input.fpkm_tracking> <output.gct>
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
    
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                output = value
            
        if len(args) < 2:
            raise Usage(help_message)
        
        fpkm_tracking = open(args[0])
        gct_out = open(args[1], "w")
        
        num_descriptor_cols = 6
        sample_names = []
        while True:
            h = fpkm_tracking.readline()
            h = h.strip()
            if h != "":
                h_cols = h.split("\t")
                if len(h_cols) < 12:
                    print >> sys.stderr, "Error: malformed header"
                    sys.exit(1)
                num_samples = (len(h_cols) - num_descriptor_cols) / 3
                for i in range(0, num_samples):
                    FPKM_label = h_cols[num_descriptor_cols + (3 * i)]
                    name_end = FPKM_label.rfind("_FPKM")
                    name = FPKM_label[:name_end]
                    sample_names.append(name)
                break 
        
        #print sample_names
        expr_records = []
        for line in fpkm_tracking.readlines():
            line = line.strip()
            if line == "":
                continue
            cols = line.split("\t")
            if len(cols) != 6 + (3 * num_samples):
                print >> sys.stderr, "Error: malformed record"
                sys.exit(1)
            expr_strs = []
            for i in range(0, num_samples):
                FPKM_string = cols[num_descriptor_cols + (3 * i)]
                FPKM = float(FPKM_string)
                if math.isnan(FPKM) or math.isinf(FPKM):
                    FPKM = 0.0
                FPKM = str(FPKM)
                expr_strs.append(FPKM)
            rec = []
            rec.append(cols[0])
            desc = "na|@%s|" % cols[5]
            rec.append(desc)
            rec.extend(expr_strs)
            expr_records.append(rec)
        
        print >> gct_out, "#1.2"
        print >> gct_out, "%d\t%d" %  (len(expr_records), num_samples)
        print >> gct_out, "NAME\tDescription\t%s" % ("\t".join(sample_names))   
        for rec in expr_records:
            print >> gct_out, "%s" % "\t".join(rec) 
        
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())
