#! /bin/bash



python meop_publish.py --copydata
python meop_publish.py --genplots
python meop_publish.py --genmaps
python meop_publish.py --compress --rebuild

