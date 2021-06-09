#!/bin/bash

function="try, $1; catch; end; quit;"
echo $function
filename=log_process_function_$(date +"%Y%m%d_%H%M%S").txt

vncserver :2
matlab -nosplash -nodesktop -logfile $filename -display :2 -r "${function}"

exit;
