#!/bin/bash

filename1='Tuttle_pdbseq.csv'
filename2='Tuttle_pdbseq_NoHyphen.csv'


awk -F, '{print $1}' Tuttle_AP_50ns.csv > Tuttle_pdbseq.csv

tr "-" " "< Tuttle_pdbseq.csv > Tuttle_pdbseq_NoHyphen.csv
