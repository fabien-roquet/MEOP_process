#!/bin/bash

# process a tag
#python meop_process.py --smru_name ft23-C899-18 --doc_latex

# process a deployment
# DEPL='ct161'
#python meop_process.py --deployment ft23 --doc_latex

# create figures
#python meop_process.py --deployment ft23 --plots --doc_latex



# TASK='--do_all'
# DEPL='ft23'
TASK='--doc_latex'
#  python meop_process.py --deployment wd3 ${TASK}
#  python meop_process.py --deployment ft23 ${TASK}
#  python meop_process.py --deployment ct81 ${TASK}
#  python meop_process.py --deployment ct125 ${TASK}
#  python meop_process.py --deployment ct104 ${TASK}
# python meop_process.py --deployment ct134 ${TASK}
# python meop_process.py --deployment ct10 ${TASK}
#  python meop_process.py --deployment ct107 ${TASK}
#  python meop_process.py --deployment ct160 ${TASK}
#  python meop_process.py --deployment ct161 ${TASK}
#  python meop_process.py --deployment ct164 ${TASK}
#  python meop_process.py --deployment ct166 ${TASK}
#  python meop_process.py --deployment ct168 ${TASK}
#  python meop_process.py --deployment ct171 ${TASK}
#  python meop_process.py --deployment ct161 ${TASK}

# python meop_process.py --deployment ct158 ${TASK}
# python meop_process.py --deployment ct155 ${TASK}
# python meop_process.py --deployment ct159 ${TASK}
# python meop_process.py --deployment ct162 ${TASK}
# python meop_process.py --deployment ct163 ${TASK}
# python meop_process.py --deployment ct165 ${TASK}
# python meop_process.py --deployment ${DEPL} ${TASK}

#12/06/23
# python meop_process.py --deployment ct175 ${TASK}
# python meop_process.py --deployment ct177 ${TASK}
# python meop_process.py --deployment ct156 ${TASK}
# python meop_process.py --deployment ct157 ${TASK}
# python meop_process.py --deployment wd2 ${TASK}
# python meop_process.py --deployment wd22 ${TASK}
# python meop_process.py --deployment ct170 ${TASK}
# python meop_process.py --deployment ft33 ${TASK}
# python meop_process.py --deployment ft34 ${TASK}
# python meop_process.py --deployment wd27 ${TASK}
# python meop_process.py --deployment ct176 ${TASK}

python meop_process.py --do_all
# python meop_postprocess.py --do_all



python meop_metadata.py --do_all

# rebuild lprofiles
# python meop_metadata.py --rebuild --qf lr0
# python meop_metadata.py --rebuild --qf lr1
# python meop_metadata.py --rebuild --qf hr0
# python meop_metadata.py --rebuild --qf hr1
# python meop_metadata.py --rebuild --qf fr0
# python meop_metadata.py --rebuild --qf fr1



python meop_publish.py --do_all

# python meop_publish.py --copydata
# python meop_publish.py --global_attributes
# python meop_publish.py --create_list_profile
# python meop_publish.py --genplots
# python meop_publish.py --genmaps
# python meop_publish.py --compress


