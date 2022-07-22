# 4DNucleome
4DNucleome theoretical analyzes and numerical simulates dynamic behaviors of transcriptional bursting across space and time described in the paper "Fundamental Principles of Enhancer-Promoter Communication in Transcriptional Bursting".
## **Environment setup**

First let’s set up the environment to run the simulations. The simulated system is complex and requires parallel computing. We can put the simulated code on the server for computing, or test a small demo on the local machine. First of all, download the repository as a zip file clicking the green button in the top right of this page to your local machine. 

## **Environment setup**

Let’s set up the environment to run the simulations. The simulated system is complex and requires parallel computing. We can put the simulated code on the server for computing, or test a small demo on the local machine. Here, we briefly introduce the steps of setting up the environment locally.

### **Download the scripts**

First of all, download the repository as a zip file clicking the green button in the top right of this page to your local machine. 

### **Requirements**

The software prerequisites for the package to work are Matlab R2018a or later version of Matlab. You can start parallel pool (parpool) using the command `parpool`. Then, you can obtain the number of workers. 

```matlab
%% current parallel pool
if isempty(gcp('nocreate'))
    parpool(4); % 4 is the number of workers.It can be adjusted according to different local devices and servers
end
```

### **Loading the scripts**

The next step is to manually add the scripts in the same folder.

## **Example: Running, analyzing and drawing a simulation demo**

You could run ‘*testDemo.m*’ which will run a simulation. The script gives a small demo of our model. The demo contains the numerical simulation, the data analysis and the theoretical analysis as well as the drawing process. The data will generate a .mat file and store in the folder generated according to the parameters.

## Directories

The folder contains several .m files. Next, we give detailed explanation each script:

#### testKEP.m, testdG.m

The script can simulate and analyze different values of E-P communication strength $k_{EP}$ (or E-P genomic distance $d_G$). For a certain $k_{EP}$ (or $d_G$), we simulate multiple times in parallel and store the simulated data in the folder generated according to the parameters.

#### ParametersBurst.m

This function generates a parameter structure to be passed to the simulation framework. The parameters are saved fields of a structure called `params`. We need to pass in the user defined parameters and set some default parameter values.

#### SimulateBurst.m

This .m files is core simulation framework. This script includes pre-allocating variable, initializing model, updating model and saving results. 

#### InitializeConnectivityMatrix.m

This function is used to initialize generalized Rouse model, and was designed for initialized params structure input format. We supply this file to generate connectivity matrix to represent the connection of linear monomers. This file should not be used for other purposes. 

#### AnalyseBurst.m

This .m file is the core analysis process. It contains two parts: Statistical analysis of simulation data and theoretical calculation of statistical indicators. We can obtain the results numerically and theoretically containing: burst size, on-state dwell time, off-state dwell time, cycle time, burst frequency, etc.

#### AnalyseBurstPDF.m

This .m file can theoretically calculate the probability density function of burst size, on-state dwell time, off-state dwell time and cycle time.

#### MutualInformationKEP.m, MutualInformationdG.m

The script can calculate mutual information of E-P communication strength $k_{EP}$ (or  $d_G$) and transcriptional burst rates ($k_{on1}$, $k_{rec}$, $k_{rel}$). 

#### drawApproximationKEP.m, drawApproximationdG.m

This script uses a deterministic binary rate for transcriptional burst to approximately calculate burst size and burst frequency. 

#### drawPartialDerivativeKEP.m, drawPartialDerivativedG.m

This script theoretically computes the logarithmic gains (i.e., two partial derivatives) of burst size and burst frequency with respect to $k_{EP}$ (or $d_G$).

#### drawApproximationSlopeHeatmap.m, drawApproximationSlopeSperation.m

This script theoretically maps a high dimensional parameter space into an experimentally measurable and theoretically computable two-dimensional space, and then calculate the ratio of $\frac{S_{BS}}{S_{BF}}$  to show which of burst size and burst frequency is primarily regulated by E-P communication. Finally, we show the region that the E-P communication mainly modulates burst frequency rather than burst size. 

#### drawParamEffectKEP.m, drawParamEffectdG.m, drawParamPowerlaw.m

This script theoretically investigates the effects of model parameters on burst size and burst frequency with respect to $k_{EP}$ (or $d_G$) and shows the power-law behaviors.

#### EPSpatDistPDFKEP.m, EPSpatDistPDFdG.m, drawChromatinDyns.m

`EPSpatDistPDFKEP.m` (or `EPSpatDistPDFdG.m`) simulated the E-P spatial dynamics and theoretically compute the distribution of E-P spatial distance. `drawChromatinDyns.m` plots curves of upstream chromatin dynamics.

#### drawKEP.m, drawdG.m, drawPDF.m, drawCDF.m, drawEnhancerDelete.m, drawTimescaleSperation.m

These .m files are used to draw figure. `drawKEP.m` (or `drawdG.m`) plots curves of the effect of  $k_{EP}$ (or $d_G$) on burst size and burst frequency. `drawPDF.m` draws the probability density function of burst size, on-state dwell time, off-state dwell time, cycle time. `drawCDF.m` draws cumulative distribution function of different $k_{EP}$ or $d_G$. `drawEnhancerDelete.m` draws the effect of enhancer deletion. `drawTimescaleSperation.m` draws the effect of model parameters to scale speration.

#### drawTRKEP.m, drawTRdG.m

This script draws the traveling ratio of  $k_{EP}$ (or $d_G$) .

