# 2D_ALDIC 
AL-DIC(Augmented Lagrangian DIC) is a fast, parallel-computing hybrid DIC algorithm, which combines advantages of local subset DIC method (fast computation speed, and parallel computing) and finite-element-based global DIC method (guarantee global kinematic compatibility and decrease noise).  

## Advantages of AL-DIC algorithm:
* [1] It’s a fast algorithm using distributed parallel computing.  
* [2]	Global kinematic compatibility is added as a global constraint in the form of augmented Lagrangian, and solved using Alternating Direction Method of Multipliers scheme.
* [3]	Both displacement fields and affine deformation gradients are correlated at the same time.
* [4]	No need of much manual experience about choosing displacement smoothing filters.
* [5]	It works well with compressed DIC images and adaptive mesh. See our paper: Yang, J. & Bhattacharya, K. Exp Mech (2019). https://doi.org/10.1007/s11340-018-00459-y;
* [6]	Both accumulative and incremental DIC modes are implemented to deal with image sequences, which is especially quite useful for very large deformations.

## Usage and citation:
* [1] Yang, J. (2019, March 6). 2D_ALDIC (Version 1.0). CaltechDATA. https://doi.org/10.22002/d1.1188
* [2] Yang, J. and Bhattacharya, K. Augmented Lagrangian Digital Image Correlation. Exp.Mech. 59: 187, 2018. https://doi.org/10.1007/s11340-018-00457-0.
 


