# dscr-template
a template repository for a dynamic statistical comparison

# How to use this to set up a new DSC

1. Clone this repo.
2. Put at least one datamaker function in a .R file in the datamakers directory (all .R files in this directory will be sourced before scenarios.R). See the file datamakers/eg_datamaker.R for example.
3. Put at least one method function in a .R file in the methods directory (all .R files in this directory will be sourced before methods.R). See the file methods/eg_method.R
4. edit the file scenarios.R to define your scenarios 
5. edit the file methods.R to define your methods
6. edit the file score.R to define your scoring function
7. replace the text in this README.md file with a description of the DSC. Include background, and definitions of the structure of the objects 'meta', 'input', and 'output' that is used by your DSC.
8. Run your DSC by running 'source("run_dsc.R")' in R. [Make sure you have installed the 'dscr' package first!]
