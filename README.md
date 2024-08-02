# SBS-GEnII-Analysis
Collection of scripts for SBS GEn-II analysis

## To use the python scripts:
Add the following line to your login script:
```
module load python/2.7.18
```
Go to your SBS-GEnII-Analysis directory and type the following:
```
source setup_py.csh
```
This will compile a python virtual environment. After you do this once, ever time you want to execute python scripts, you can do the following:
```
source sourceme_py.csh
```
Sourcing either of these files will put you in a virtual environment. You can tell because the start of your terminal command lines will begin with "(venv)". You can code as normal within this, but to exit this environment, simply type the following:
```
deactivate
```
