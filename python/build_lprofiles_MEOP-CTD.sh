#! /bin/bash


# rebuild lprofiles
python meop_metadata.py --rebuild --qf lr0
python meop_metadata.py --rebuild --qf lr1
python meop_metadata.py --rebuild --qf hr0
python meop_metadata.py --rebuild --qf hr1
python meop_metadata.py --rebuild --qf fr0
python meop_metadata.py --rebuild --qf fr1

