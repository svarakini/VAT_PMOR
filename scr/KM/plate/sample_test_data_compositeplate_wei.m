addpath(genpath([pwd filesep 'scr' filesep 'FEM_pre']));

addpath(genpath([pwd filesep 'scr' filesep 'KM']));
%
addpath(genpath([pwd filesep 'scr' filesep 'Boundary_conditions']));
addpath(genpath([pwd filesep 'scr' filesep 'FORCE']));
addpath(genpath([pwd filesep 'scr' filesep 'mis']));
addpath(genpath([pwd filesep 'scr' filesep 'FEM_post']));
%
addpath(genpath([pwd filesep 'scr' filesep 'response']));
addpath(genpath([pwd filesep 'scr' filesep 'solver']));
addpath(genpath([pwd filesep 'scr' filesep 'VAT_para']));

addpath(genpath([pwd filesep 'scr' filesep 'export_fig']));

addpath(genpath([pwd filesep 'scr' filesep 'mesh']));



% MatProperty;
psi2pa=6894.75729;
% psipa=1;
Mat.kappa=5/6;
Mat.E1=181e9;
Mat.E2=10.27e9;
Mat.G12=7.17e9;
Mat.G13=4e9;
Mat.G23=4e9;
Mat.v12=0.28;
Mat.v21=Mat.E2/Mat.E1*Mat.v12;
Mat.density=1800;
Mat.alpha1=-0.04e-6;
Mat.alpha2=16.7e-6;

% to from a 3D orthotropic material
Mat.E3  = Mat.E1;
Mat.v13 = Mat.v12;
Mat.v31 = Mat.v13/Mat.E1*Mat.E3;
Mat.v23 = Mat.v12;
Mat.v32 = Mat.v23*Mat.E3/Mat.E2;

%  T01 = [-30 60];
%% Laminate configuration
half_number = 4;
flag = 'SYM'; % symmetric laminate
 
Laminate.layer_thickness =1.272e-4;
Laminate.layer= half_number*2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------PLATE STIFFNESS Matrix --------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Laminate.thickness=Laminate.layer_thickness*ones(1,Laminate.layer);%% Uniform thickness

%% Total Material Stiffness and Gometric Stiffness, Mass matrix;

Plate.length = .30;
Plate.width  = .30;

post_plot = 0; % 0 - disable plot; 1 - enable plot option

%%%%%% % Stru.bdfname ='generate_mesh_8noded.bdf';
%% NASTRUAL TO PHYSICAL

FEM = mesh_QUAD8_v2(8,8); % mesh a 2D natural surface [-1,1]

FEM.numberElements = size(FEM.elementNodes ,1);
FEM.typeplate='CQUAD8';
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

%% =============== PANEL GEOMETRY================================
%%%==============================================================
FEM.nodesCord(:,2) = FEM.nodesCord(:,2)/max(FEM.nodesCord(:,2))*Plate.length ;
FEM.nodesCord(:,3) = FEM.nodesCord(:,3)/max(FEM.nodesCord(:,3))*Plate.width ;
%---------Coordinates of node------------
Xcoord = FEM.nodesCord(:,2);
Ycoord = FEM.nodesCord(:,3);
Zcoord = FEM.nodesCord(:,4);

Stru.length = max(abs(Xcoord));
Stru.width  = max(abs(Ycoord));

%
Stru.thickness= sum(Laminate.thickness);

% The fiber path orientation is generalized for curvilinear fiber
% determine the middle point in fiber path angle computation
center =[Stru.length Stru.width]/2;

% length of the plate
physical_length = Stru.length; % change along width direction


%% ============= External FORCE ==================
FEM.Fz=0;% p is the uniform pressure.
FEM.Mx=0;
FEM.My=0;
%%%
FEM.p= -FEM.Fz;

%% ============ IN-PLANE LOADS (N) ==========
Alpha= 0; % Alpha - load type
RHS_nodes = find(FEM.nodesCord(:,2) == Plate.length);

for iii = 1:length(RHS_nodes)
    y = FEM.nodesCord(RHS_nodes(iii),3);
    force_temp(RHS_nodes(iii),1) =  (1-Alpha*y/Plate.width) ;
end

N0= 1e3;%%

if Alpha== 2
    force_function = @(width_y) N0*(1-Alpha*width_y/Plate.width);
    total_force_applied =  integral(force_function,0,Plate.width/2);
    
    force_temp = -force_temp/sum(force_temp(RHS_nodes(1:(length(RHS_nodes)-1)/2,1)))* total_force_applied; %%25.5
else
    bot_load = 1; top_load = 1-Alpha;

    force_function = @(width_y) N0*(1-Alpha*width_y/Plate.width);
    total_force_applied =  integral(force_function,0,Plate.width);

    force_total = sum(linspace(top_load,bot_load,101)); % 101 is the nastran model nodes in right edge
    
    force_temp = -force_temp/sum(force_temp) *  total_force_applied ;
end

T01 =[45 0];
T01 = [-7.5 45];
T01 = [17.5 17.5];
T01 = [-7.5 50]; % Variable angle tow laminates
T01 =[27.5 27.5]; % straight fiber path

T0T1 = VAT_fiber(T01,half_number,flag);

% % % % QI case
% % % T0T1 = [45 45
% % %     -45 -45
% % %     0 0
% % %     90 90
% % %     90 90
% % %     0 0
% % %     -45 -45
% % %     45 45];