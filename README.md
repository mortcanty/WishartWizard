The Wishart Wizard is a graphical user interface (GUI) implemented in IDL as an extension
for the remote sensing image analysis environment ENVI. It provides a simplified and 
user-friendly platform for performing multivariate change detection with bitemporal polarimetric SAR imagery.
The change detection procedure implemented exploits the complex Wishart distribution of 
polarimetric SAR image observations in look-averaged covariance matrix format in order 
to define a per-pixel change/no-change hypothesis test. It includes 
approximations for the probability distribution of the
test statistic, and so permits quantitative significance levels
to be quoted for change pixels. In addition, an improved multivariate method
is used to estimate the equivalent number of looks (ENL) of the look-averaged images,
which is a critical parameter of the hypothesis test.

If SARscape is licensed in the ENVI environment, the Wizard accesses the SARscape API at the
IDL scripting level to expose only that functionality necessary for the change detection 
analysis, thus guiding the user and simplifying the  processing chain. If SARscape is not 
present, open source alternatives for the required preprocessing are available and are 
described in this documentation.
