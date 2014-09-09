############################################################################################################
# DESCRIPTION:                        This script creates ashr.no.cxx.tar.gz, a version of the ashr package 
#                                     that does not require Rcpp
# USAGE:                              sh create.ashr.no.cxx.sh
# NOTE:                               Folder ashr (containing the R, src, ... folders) should be in same 
#                                     directory as script create.ashr.no.cxx.sh
############################################################################################################

#make a copy of ashr in ashr.no.cxx
mkdir ashr.no.cxx
cp -r ashr ashr.no.cxx/ashr

#replace cxx=TRUE with cxx=FALSE in file ash.R
cat ashr.no.cxx/ashr/R/ash.R | sed 's/cxx[ ]*=[ ]*TRUE/cxx=FALSE/' > ashr.no.cxx/ashr/R/tmp.ash.R
mv ashr.no.cxx/ashr/R/tmp.ash.R ashr.no.cxx/ashr/R/ash.R

#remove src folder
rm -r ashr.no.cxx/ashr/src

#remove lines that contain Rcpp in file DESCRIPTION 
#note: this will remove *every* line in file DESCRIPTION that contain the string "Rcpp"
cat ashr.no.cxx/ashr/DESCRIPTION | grep -v Rcpp > ashr.no.cxx/ashr/tmp.DESCRIPTION
mv ashr.no.cxx/ashr/tmp.DESCRIPTION ashr.no.cxx/ashr/DESCRIPTION

#remove lines that contain useDynLib in file NAMESPACE
cat ashr.no.cxx/ashr/NAMESPACE | grep -v useDynLib > ashr.no.cxx/ashr/tmp.NAMESPACE
mv ashr.no.cxx/ashr/tmp.NAMESPACE ashr.no.cxx/ashr/NAMESPACE

#remove RcppExports.R file
rm ashr.no.cxx/ashr/R/RcppExports.R

#create tar archive
cd ashr.no.cxx
tar -pczf ashr.no.cxx.tar.gz ashr/
mv ashr.no.cxx.tar.gz .. 
cd .. 

#clean
rm -r ashr.no.cxx
