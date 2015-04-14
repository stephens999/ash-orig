# make figures
# make dsc
 
all: $(TEXFILE).pdf $(OUT_FILES) $(CROP_FILES)

figures: figures/egcdf.pdf figures/mean_cdf.pdf figures/pi0_est.pdf figures/scenarios_density.pdf
 
dsc: ../dsc-shrink/res.RData ../dsc-stability/res.RData ../dsc-znull/res.RData

../dsc-shrink/res.RData:
	cd ../dsc-shrink; R CMD BATCH run_dsc.R
 
#figures
figures/egcdf.pdf: ../dsc-shrink/res.RData
	cd ../dsc-shrink; R CMD BATCH plot_cdf_eg.R
figures/mean_cdf.pdf: ../dsc-shrink/res.RData
	cd ../dsc-shrink; R CMD BATCH plot_cdf_eg.R
figures/scenarios_density.pdf: ../dsc-shrink/res.RData
	cd ../dsc-shrink; R CMD BATCH plot_egdens.R
figures/pi0_est.pdf: ../dsc-shrink/res.RData
	cd ../dsc-shrink; R CMD BATCH plot_pi0_est.R

 
 
# clean up all files; use with caution
superclean: 
	cd ../dsc-shrink; rm -r dsc-shrink-files
	cd ../dsc-stability; rm -r dsc-stability-files

	rm figures/*.pdf

 
.PHONY: all figures superclean
