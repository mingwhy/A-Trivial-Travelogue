nicepath(https://github.com/EPFL-LCSB/nicepath) used 'git lfs', need to install git lfs first.
$ brew install git-lfs

then fetch nicepath via git command line to correctlyl download all git lfs files.
$ git clone https://github.com/EPFL-LCSB/nicepath.git

$ cd /path/to/nicepath
$ make

need to create 'output/' folder first before running below.
$ cd nicepath
$ python main.py KEGG 
# input /Users/mingyang/Downloads/nicepath/input/KEGG
# output /Users/mingyang/Downloads/nicepath/output/KEGG

$ python main.py Test 
# input /Users/mingyang/Downloads/nicepath/input/Test
# output /Users/mingyang/Downloads/nicepath/output/Retrobio

look at the `parameters.txt` file for different settings.

modify:
output_folder_path|../../nicepath/output/KEGG/
target_compounds|C01197
# List of precursor compounds: can be an entry number, a list of entry numbers as a .txt file, or 'model' to search against all model compounds
precursor_compounds|C00082


######################
conda create -n "network" python=3.8 ipython networkx==2.3
conda activate network
#conda deactivate 

######################
`ATLASxAnalyses` from Expanding biochemical knowledge and illuminating metabolic dark matter with ATLASx
https://github.com/EPFL-LCSB/ATLASxAnalyses
Download repository
$ git clone https://github.com/EPFL-LCSB/ATLASxAnalyses

To install the required dependencies:
$ cd ATLASxAnalysis $ make
# if `make` not working, manually install 
$ pip install matplotlib
$ pip install seaborn
$ pip install pandas

Reproduce Network Analysis
$ cd NetworkAnalysis/Source

Plot component distribution for database scopes
$ python3 get_component_distribution.py #Runtime: 97s

pip install --force-reinstall numpy==1.20.3


################################################################ 
the reason for `python=3.8 ipython networkx==2.3` and `numpy==1.20.3` 

ImportError: this version of pandas is incompatible with numpy < 1.20.3
your numpy version is 1.20.0.
Please upgrade numpy to >= 1.20.3 to use this pandas version

debug: AttributeError: module 'networkx' has no attribute 'read_gpickle'
solution: https://github.com/rkistner/chinese-postman/issues/21
`it's still present up till 2.3, and removed in 2.4.`
I've tried
$ pip install --force-reinstall networkx==2.3
but this induced other problems.

$ pip install --force-reinstall networkx
Installing collected packages: networkx
  Attempting uninstall: networkx
    Found existing installation: networkx 2.3
    Uninstalling networkx-2.3:
      Successfully uninstalled networkx-2.3
Successfully installed networkx-3.1



(network) Source $pip list
Package             Version
------------------- ------------
appnope             0.1.2
asttokens           2.0.5
backcall            0.2.0
contourpy           1.1.1
cycler              0.12.1
decorator           5.1.1
executing           0.8.3
fonttools           4.46.0
importlib-resources 6.1.1
ipython             8.12.2
jedi                0.18.1
kiwisolver          1.4.5
matplotlib          3.7.4
matplotlib-inline   0.1.6
networkx            2.3
numpy               1.20.3
packaging           23.2
pandas              2.0.3
parso               0.8.3
pexpect             4.8.0
pickleshare         0.7.5
Pillow              10.1.0
pip                 23.3.1
prompt-toolkit      3.0.36
ptyprocess          0.7.0
pure-eval           0.2.2
Pygments            2.15.1
pyparsing           3.1.1
python-dateutil     2.8.2
pytz                2023.3.post1
seaborn             0.13.0
setuptools          68.0.0
six                 1.16.0
stack-data          0.2.0
traitlets           5.7.1
typing_extensions   4.7.1
tzdata              2023.3
wcwidth             0.2.5
wheel               0.41.2
zipp                3.17.0





