
all: modeSolver4l_LpSl_Verbose_mex.mexa64 fo4l_LpSl_mex.mexa64 field4l_LpSl_Verbose_mex_ext.mexa64

modeSolver4l_LpSl_Verbose_mex.mexa64: modeSolver4l_LpSl_Verbose_mex.cpp dispersionLP.cpp DmatrixLp.cpp SmatrixLP.cpp gen_deltan_prof.cpp gen_rho.cpp modeSolverLP.cpp refIndProfile.cpp sellmeierEquation.cpp zeroFinderF.cpp zeroFinderGE.cpp zeroFinderMode.cpp
	mex -lgsl -lgslcblas $^

fo4l_LpSl_mex.mexa64: fo4l_LpSl_mex.cpp dispersionLP.cpp DmatrixLp.cpp SmatrixLP.cpp FO_dmd.cpp gen_deltan_prof.cpp gen_rho.cpp modeSolverLP.cpp refIndProfile.cpp sellmeierEquation.cpp zeroFinderF.cpp zeroFinderGE.cpp zeroFinderMode.cpp
	mex -lgsl -lgslcblas $^

field4l_LpSl_Verbose_mex_ext.mexa64: field4l_LpSl_Verbose_mex_ext.cpp fieldLP.cpp dispersionLP.cpp DmatrixLp.cpp SmatrixLP.cpp FO_dmd.cpp gen_deltan_prof.cpp gen_rho.cpp modeSolverLP.cpp refIndProfile.cpp sellmeierEquation.cpp zeroFinderF.cpp zeroFinderGE.cpp zeroFinderMode.cpp
	mex -lgsl -lgslcblas $^

.PHONY: clean

clean:
	rm *.mexa64
