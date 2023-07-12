cd /project/st-ljfoster-1
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh
# /project/st-ljfoster-1/miniconda3

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/project/st-ljfoster-1/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/project/st-ljfoster-1/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/project/st-ljfoster-1/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/project/st-ljfoster-1/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

mkdir /project/st-ljfoster-1/nps-screening
cd !$

conda create --prefix /project/st-ljfoster-1/nps-screening/env
conda activate /project/st-ljfoster-1/nps-screening/env

conda install r-tidyverse r-magrittr r-argparse
conda install -c bioconda bioconductor-xcms
