# Multisp_model_JAE

Model and data used in the following paper:  

**Trijoulet V., Fay G. and Miller T.J. 2019 Performance of a state-space multispecies model: what are the consequences of ignoring predation and process errors in stock assessments? Journal of Applied Ecology**

#### Scripts:  
* The estimation models (EMs) can be fitted from the "**Fitting_the_EMs.R**" file  
* The structural code is in the file "**MS_SSM.cpp**"  
* Functions used in the cpp code are in "**MS_SSM_functions.h**"  
* The script "**MS_SSM_make.map.fn.R**" is the function that creates the map argument of the objective function depending on the options chosen in the models  
* The script "**MS_SSM_make.random.fn.R**" is the function that creates the random argument of the objective function depending on the options chosen in the models  

#### Data:
* The 1000 data sets used in the paper are available as lists in "**sim_data_all.RData**"  
* A set of initial values is available in "**init.RData**"  

#### Please note:
This model is not the most recent version of the multispecies model and it appears as it was used in the paper. 
A new and improved version of the model is available in the repository "MS_SSM". Please consider this latter repository for an updated version of the model.