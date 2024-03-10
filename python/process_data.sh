#!/bin/bash

# process a tag
#python meop_process.py --smru_name ft23-C899-18 --doc_latex

# process a deployment
# DEPL='ct161'
#python meop_process.py --deployment ft23 --doc_latex

# create figures
#python meop_process.py --deployment ft23 --plots --doc_latex



#  TASK='--do_all'
# DEPL='ft23'
# TASK='--do_all'
# python meop_process.py --deployment tu116 ${TASK}
# python meop_process.py --deployment tu117 ${TASK}
# python meop_process.py --deployment wd22 ${TASK}
# python meop_process.py --deployment tu120 ${TASK}
# python meop_process.py --deployment wd27 ${TASK}
# python meop_process.py --deployment wd31 ${TASK}
#  python meop_process.py --deployment wd19 ${TASK}
#  python meop_process.py --deployment wd20 ${TASK}
#  python meop_process.py --deployment wd30 ${TASK}
#  python meop_process.py --deployment wd31 ${TASK}


# python meop_process.py --deployment ct160 ${TASK}
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

# python meop_process.py --deployment ct172 ${TASK}
# python meop_process.py --deployment ct173 ${TASK}
# python meop_process.py --deployment ct177 ${TASK}
# python meop_process.py --deployment ct179 ${TASK}
# python meop_process.py --deployment ct181 ${TASK}
# python meop_process.py --deployment ct72 ${TASK}

# python meop_process.py --do_all
# python meop_postprocess.py --do_all


python meop_process_all.py 
python meop_metadata.py --rebuild

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


