
# How do I install Gurobi for Matlab?
# https://support.gurobi.com/hc/en-us/articles/4533938303505-How-do-I-install-Gurobi-for-Matlab-
the default <installdir> for Gurobi 10.0.2 is /Library/gurobi1002/macos_universal2 for Mac
To get started, type the following commands within MATLAB to change to the matlab directory and call gurobi_setup:

>> cd /Library/gurobi1002/macos_universal2/matlab
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
or you need to type below everytime you use
>> addpath('/Users/ming/Documents/bioinfo_software/GEM_metabolic.models/RAVEN-main')

view all path:
>> path
#####################################################
# install gurobi (solver algorithm)

# https://www.gurobi.com/features/academic-named-user-license/
# https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation_guide.html
# https://support.gurobi.com/hc/en-us/articles/14462206790033-How-do-I-install-Gurobi-for-R-
# https://sites.google.com/view/jessica-leung/home/teaching/using-gurobi-with-r-and-python
##Using Gurobi with RStudio
#1. Register a Gurobi account as an academic user at https://pages.gurobi.com/registration . Please register as a student using your university email address and don’t lose the password.
#2. Go to https://www.gurobi.com/downloads/gurobi-optimizer-eula/ , read and accept the End User License Agreement. Then download the current version of Gurobi optimizer. 
# Install the optimizer with default settings.
#3. Request for a license at  https://www.gurobi.com/downloads/end-user-license-agreement-academic/, enter User Portal, and you will obtain a Gurobi key. 
# Copy and paste the Gurobi key to the terminal: $grbgetkey xxxxxx
#4. Run the following R command in RStudio:
install.packages('/Library/gurobi1002/macos_universal2/R/gurobi_10.0-2_R_4.2.0.tgz', repos=NULL)
#5. Now you should be able to solve linear programming problem with Gurobi in RStudio. Try running the following piece of code:
  
# This example formulates and solves the following simple LP model:
#  maximize
#        x + 2 y + 3 z
#  subject to
#        x +   y <= 1
#        y +   z <= 1
library(Matrix)
library(gurobi)
M <- list()
M$A<- matrix(c(1,1,0,0,1,1), nrow=2, byrow=T)
M$obj<- c(1,2,3)
M$modelsense <- 'max'
M$rhs <- c(1,1)
M$sense <- c('<', '<')
result <- gurobi(M)
print(result$objval)
print(result$x)
