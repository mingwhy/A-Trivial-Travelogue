transcript2flux
===============

Scripts and results from gene expression methods assessment paper:

D. Machado, M.J. Herrgard. Systematic Evaluation of Methods for Integration of Transcriptomic Data into Constraint-Based Models of Metabolism. PLoS Computational Biology, 10(4), 2014.

http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003580


INSTRUCTIONS:

- Install cobra toolbox for matlab from: http://opencobra.sourceforge.net/
- Install MADE from: http://www.bme.virginia.edu/csbl/downloads-made.php
- Install a solver with support for LP, QP and MILP (recommended: Gurobi 5.5)

- Use main.m to run all tests from the paper.

# install cobra toolbox from https://opencobra.github.io/cobratoolbox/stable/installation.html

cd  /Users/ming/Documents/bioinfo_software/GEM_metabolic.models 
git clone --depth=1 https://github.com/opencobra/cobratoolbox.git cobratoolbox
#Change to the folder cobratoolbox/ and run from matlab
>> initCobraToolbox

# set up matlab interface to gurobi
#https://www.gurobi.com/documentation/current/quickstart_mac/matlab_interface.html
#https://www.gurobi.com/documentation/current/quickstart_mac/matlab_setting_up_grb_for_.html
 To get started, type the following commands within MATLAB to change to the matlab directory and call gurobi_setup:

>> cd /Library/gurobi1002/macos_universal2/matlab
>> gurobi_setup
>> gurobi_setup

The MATLAB interface for Gurobi 10.0.2 has been installed.

The directory
    /Library/gurobi1002/macos_universal2/matlab/
has been added to the MATLAB path.
To use Gurobi regularly, you must save this new path definition.
To do this, type the command
    savepath
at the MATLAB prompt. Please consult the MATLAB documentation
if necessary.
>> savepath

#run example: https://www.gurobi.com/documentation/current/quickstart_mac/matlab_example.html
>> cd /Library/gurobi1002/macos_universal2/examples/matlab
>> mip1




