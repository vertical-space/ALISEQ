#!/bin/bash
# properties = {properties}
 
source $HOME/.bashrc
 
export PATH=$PATH:/exports/cmvm/eddie/eb/groups/watson_grp/software/
 
module add anaconda 1>/dev/null 2>/dev/null
 
{exec_job}

