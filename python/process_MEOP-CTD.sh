#! /bin/bash


# process a tag
#python meop_process.py --smru_code ft23-C899-18

# process a deployment
#python meop_process.py --deployment ft23

# create figures
#python meop_process.py --deployment ft23 --plots




#TASK='--do_all'
#TASK='--create_hr2'
#python meop_process.py --deployment wd3 ${TASK}
#python meop_process.py --deployment ft23 ${TASK}
#python meop_process.py --deployment ct81 ${TASK}
#python meop_process.py --deployment ct125 ${TASK}
#python meop_process.py --deployment ct104 ${TASK}
#python meop_process.py --deployment ct134 ${TASK}
#python meop_process.py --deployment ct10 ${TASK}
#python meop_process.py --deployment ct107 ${TASK}

#python meop_postprocess.py --do_all

python meop_publish.py --do_all

#python meop_publish.py --copydata
#python meop_publish.py --global_attributes
#python meop_publish.py --create_list_profile
#python meop_publish.py --genplots
#python meop_publish.py --genmaps
#python meop_publish.py --compress


