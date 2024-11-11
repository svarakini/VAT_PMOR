# VAT_PMOR
Parametric model order reduction of Variable Angle Tow composite plate

Reduced Order Finite Element Analysis of Composite Laminated Plate
Author: Varakini Sanmugadas (svarakini@vt.edu)
<https://github.com/svarakini>

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal with the Software without restriction, including without 
limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons 
to whom the Software is furnished to do so.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
THE USE OR OTHER DEALINGS WITH THE SOFTWARE. 
--------------------------------------------------------------------------

Full order buckling analysis is adapted from the repository below and
modified 
<https://github.com/zhaowei0566/SPAD>

- To run the buckling analysis
     Main file: a_main_Buckling_Analysis.m
  
- The variables to change:

%% Choose FOM/ROM model
FOM_ROM_flag = 'kROM_gFOM'; %'ROM'; %'FOM'; %  
% 'kROM_gFOM' - Reduced order model using affine decompositiion for linear
%               stiffness matrix
% 'ROM' - Reduced order model using affine decompositiion for both linear 
%         and geometric stiffness matix
% 'FOM' - Full order model by Wei Zhao, AIAAJ 2019
% https://github.com/zhaowei0566/SPAD

% VAT angles 
% Rows correspond to layers. Add/remove rows to change number of layers.
% Column 1 - Theta_0; Column 2 - Theta_1;
T0T1 = [45 45
    -45 -45
    0 0
    90 90
    90 90
    0 0
    -45 -45
    45 45
    ];

num_layers = size(T0T1,1);

% layer thickness
layer_thick=0.0001*ones(1,num_layers);

% load type
DisplacemntBoundaryCase='uniform';  

mat_files_folder=[pwd filesep 'mat_files_latin_vat' filesep 'case_n_train400_buck24_stat42'];

'..\..\configurations'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

