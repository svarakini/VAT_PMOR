%% Reduced Order Finite Element Analysis of Composite Laminated Plate
% Author: Varakini Sanmugadas (svarakini@vt.edu)
% <https://github.com/svarakini>
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal with the Software without restriction, including without 
% limitation the rights to use, copy, modify, merge, publish, distribute, 
% sublicense, and/or sell copies of the Software, and to permit persons 
% to whom the Software is furnished to do so.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
% OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
% THE USE OR OTHER DEALINGS WITH THE SOFTWARE. 
%--------------------------------------------------------------------------
%
% Full order buckling analysis is adapted from the repository below and
% modified 
% <https://github.com/zhaowei0566/SPAD>
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

global Mat FEM  Stru Laminate Plate Laminates

% ==== add subroutines ======
addpath(genpath([pwd filesep 'scr' filesep 'FEM_pre']));
addpath(genpath([pwd filesep 'scr' filesep 'KM']));
addpath(genpath([pwd filesep 'scr' filesep 'Boundary_conditions']));
addpath(genpath([pwd filesep 'scr' filesep 'FORCE']));
addpath(genpath([pwd filesep 'scr' filesep 'mis']));
addpath(genpath([pwd filesep 'scr' filesep 'FEM_post']));
addpath(genpath([pwd filesep 'scr' filesep 'response']));
addpath(genpath([pwd filesep 'scr' filesep 'solver']));
addpath(genpath([pwd filesep 'scr' filesep 'VAT_para']));
addpath(genpath([pwd filesep 'scr' filesep 'export_fig']));
addpath(genpath([pwd filesep 'scr' filesep 'mesh']));
addpath(genpath([pwd filesep 'mat_files']));

post_plot = 0; % 0 - disable plot; 1 - enable plot option


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% MESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEM = mesh_QUAD8_v2(36,36); % mesh a 2D natural surface [-1,1]
FEM.numberElements = size(FEM.elementNodes ,1);

FEM.nodesCord(:,2) = FEM.nodesCord(:,2)/max(FEM.nodesCord(:,2))*Plate.length ;
FEM.nodesCord(:,3) = FEM.nodesCord(:,3)/max(FEM.nodesCord(:,3))*Plate.width ;

%---------Coordinates of node------------
Xcoord = FEM.nodesCord(:,2);
Ycoord = FEM.nodesCord(:,3);
Zcoord = FEM.nodesCord(:,4);

%-----------D.O.F for each Node--------------
FEM.PlateNodeDof = 5;%% for plate
FEM.nodeCoordinates = FEM.nodesCord(:,2:3);

FEM.GDof = FEM.PlateNodeDof*size(FEM.nodeCoordinates,1);
FEM.typeplate='CQUAD8';
FEM.Dimension='2D';

switch FEM.typeplate
    case 'CQUAD4'
        FEM.GaussPointShear='1by1';
        FEM.GaussPointBend='2by2';
    case 'CQUAD8'
        FEM.GaussPointShear='2by2';
        FEM.GaussPointBend='3by3';
        FEM.n_GaussPoints = 8;
end

%---------------------------------------
FEM.NodeNumber=size(FEM.nodesCord,1);
FEM.elementNumber=size(FEM.elementNodes,1);
FEM.nodeCoordinates_label = zeros(size(FEM.nodeCoordinates,1),4);

FEM.nodeCoordinates_label(:,2:3) = FEM.nodeCoordinates;
if post_plot~=0
    HH = figure(200);hold on;
    plot(FEM.nodeCoordinates(:,1),FEM.nodeCoordinates(:,2),'ko')
end
FEM.nodeCoordinates_label(:,1)  = 1:size(FEM.nodeCoordinates,1);

% plot mesh of the plate
if post_plot~=0
    patch_plot(FEM.elementNodes,FEM.nodeCoordinates_label,200,'skin');axis image;
end

%% %%%%%%%%%%%%%%%%%%%%%%%% PANEL GEOMETRY %%%%%%%%%%%%%%%%%%%%
Stru.length = max(abs(Xcoord));
Stru.width  = max(abs(Ycoord));
Stru.thickness= sum(layer_thick);

Stru.center =[Stru.length Stru.width]/2;

%% %%%%%%%%%%%%%%%%%%%%%% MATERIAL PROPERTY %%%%%%%%%%%%%%%%%%%
% MatProperty;
psi2pa=6894.75729;

Mat.density=1800; Mat.kappa=5/6;
Mat.E1=181e9; Mat.E2=10.27e9;
Mat.G12=7.17e9; Mat.G13=4e9; Mat.G23=4e9;
Mat.v12=0.28; Mat.v21=Mat.E2/Mat.E1*Mat.v12;
Mat.alpha1=-0.04e-6; Mat.alpha2=16.7e-6;

% to form a 3D orthotropic material
Mat.E3  = Mat.E1; Mat.v13 = Mat.v12; Mat.v31 = Mat.v13/Mat.E1*Mat.E3; 
Mat.v23 = Mat.v12; Mat.v32 = Mat.v23*Mat.E3/Mat.E2;

% Laminate.layer_thickness =layer_thick;
half_number = Laminates.plynumber/2;
Laminate.layer= Laminates.plynumber;
Laminate.thickness=layer_thick;

% T0T1 = VAT_fiber(T01,half_number,flag);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PLATE STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FEM.typeplate='CQUAD8';

% length of the plate
physical_length = Stru.length; 

FEM.DisplacemntBoundaryCase = DisplacemntBoundaryCase;

FEM.BCtype = 'SSSS-3'; %% Four simple supported sides
[FEM.ActiveDof_ssss,FEM.SideDof_ssss] = EssentialBCPlate5Dof(FEM);


%% ============= Matrix Assembly =============

% Assemble Linear Stiffness Matrix
Kplate = kstiff_assemble_rom(FOM_ROM_flag,Stru,Laminate,Mat,FEM,T0T1,mat_files_folder);
if strcmp(FOM_ROM_flag,'FOM')
    Kplate_actv = Kplate(FEM.ActiveDof_ssss,FEM.ActiveDof_ssss); 
elseif strcmp(FOM_ROM_flag,'ROM') || strcmp(FOM_ROM_flag,'kROM_gFOM')
    Kplate_actv = Kplate;
end

% Assemble Geometric Stiffness Matrix
[KGplate_actv,FEM] = gstiff_assemble_rom(FOM_ROM_flag,Stru,Laminate,Mat,FEM,T0T1,Kplate,mat_files_folder);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% BUCKLING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% EIGENVALUE ANALYSIS %%%%%%%%%%%%%%%%%%%
% tic
[V,D] = eigs(Kplate_actv,KGplate_actv,5,'sm');
% toc
[DD,modeNo]=sort(diag(D));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD FACTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------minimal load factor---------------
mineiglabel=find(min(abs(DD))==-DD);
if isempty(mineiglabel)==1
    mineiglabel=find(min(abs(DD))==DD);
end
plotmodenNo=mineiglabel;

critical_loadfactor = DD(mineiglabel);

elements_of_interest = [];
num_interest = 1;

for elem = 1:FEM.elementNumber
    nodes_4_element = FEM.elementNodes(elem,:);
    element_yes = 0;
    for node_num = 1:length(nodes_4_element)
        cord_temp =  label2cord(nodes_4_element(node_num),FEM.nodesCord);
        if abs((cord_temp(2) - Stru.length)/Stru.length)<1e-5  % x=a, RHS
            element_yes = 1;
        end
    end
    if  element_yes ==1
        elements_of_interest(num_interest) = elem;
        num_interest = num_interest+1;
    end
end

% Calculate NXX
center=Stru.center;
calculate_stress_mat_cord_noplot

%%
% Critical buckling load
Pcr = critical_loadfactor*sum(NXX(elements_of_interest).*Stru.length/length(elements_of_interest)) ;
% Critical buckling parameter
Kcr = Pcr*Stru.length^2/(Mat.E1*Stru.width*Stru.thickness^3);


