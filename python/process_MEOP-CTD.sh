#! /bin/bash


# process a tag
#python meop_process.py --smru_code ft23-C899-18

# process a deployment
#python meop_process.py --deployment ft23

# create figures
#python meop_process.py --deployment ft23 --plots


#python meop_process.py --deployment wd3 --do_all
#python meop_process.py --deployment ft23 --do_all
#python meop_process.py --deployment ct81 --do_all
#python meop_process.py --deployment ct125 --do_all

#problem with location!?
python meop_process.py --deployment ct104 --do_all

#python meop_process.py --deployment ct134 --do_all
#python meop_process.py --deployment ct10 --do_all
#python meop_process.py --deployment ct107 --do_all


python meop_publish.py --copydata --rebuild
python meop_publish.py --global_attributes
python meop_publish.py --genplots --rebuild
python meop_publish.py --genmaps --rebuild
python meop_publish.py --compress --rebuild


