% function V_mat=lamination_parameters(stru_thickness,thickness,T0T1,Element_center_X,center,physical_length)
% svarakini@vt.edu
function [V_mat,Q_mat]=lamination_parameters_elewise_array(stru_thickness,thickness,T0T1,Element_center_X,center,physical_length,Mat)
  % This function is for STRAIGHT plies (not Variable Anlgle plies. VAT must be assembled element by element)

%%
    T0 = T0T1(:,1);
    T1 = T0T1(:,2);
    %     theta = VAT_fiber_ply_angle_1D_elewise(T0,T1,Element_center_X,center,physical_length);
    theta = VAT_fiber_ply_angle_1D_elewise(T0,T1,Element_center_X,center,physical_length);
    %size(theta)
%%
    zk1=-stru_thickness/2 + thickness'.*(1:size(T0T1,1))'; %% CHECK: when LAYER THICK DIFFERENT!!!!!!!!!!!!!!
    zk =-stru_thickness/2 + thickness'.*((1:size(T0T1,1))'-ones(size(T0T1,1),1));

%%
    V0A=(zk1-zk).*ones(8,length(theta));
    V1A=(zk1-zk).*cosd(2*theta);
    V2A=(zk1-zk).*sind(2*theta);
    V3A=(zk1-zk).*cosd(4*theta);
    V4A=(zk1-zk).*sind(4*theta);

    V0B=(zk1.^2-zk.^2)./2.*ones(8,length(theta));
    V1B=(zk1.^2-zk.^2)./2.*cosd(2*theta);
    V2B=(zk1.^2-zk.^2)./2.*sind(2*theta);
    V3B=(zk1.^2-zk.^2)./2.*cosd(4*theta);
    V4B=(zk1.^2-zk.^2)./2.*sind(4*theta);

    V0D=(zk1.^3-zk.^3)./3.*ones(8,length(theta));
    V1D=(zk1.^3-zk.^3)./3.*cosd(2*theta);
    V2D=(zk1.^3-zk.^3)./3.*sind(2*theta);
    V3D=(zk1.^3-zk.^3)./3.*cosd(4*theta);
    V4D=(zk1.^3-zk.^3)./3.*sind(4*theta);
    
    %%
    V_mat.V0A = V0A;
    V_mat.V1A = V1A;
    V_mat.V2A = V2A;
    V_mat.V3A = V3A;
    V_mat.V4A = V4A;
    
    V_mat.V0B = V0B;
    V_mat.V1B = V1B;
    V_mat.V2B = V2B;
    V_mat.V3B = V3B;
    V_mat.V4B = V4B;
    
    V_mat.V0D = V0D;
    V_mat.V1D = V1D;
    V_mat.V2D = V2D;
    V_mat.V3D = V3D;
    V_mat.V4D = V4D;


    Q11=Mat.E1/(1-Mat.v12*Mat.v21);
    Q12=Mat.v12*Mat.E2/(1-Mat.v12*Mat.v21);
    Q22=Mat.E2/(1-Mat.v12*Mat.v21);
    Q66=Mat.G12;
    % Q44=Mat.G23;
    % Q55=Mat.G13;

    % NodeIndices=FEM.elementNodes; %(elem,:);
    % ele_coords=reshape(FEM.nodeCoordinates_label( NodeIndices,2),FEM.elementNumber,8);
    % centerX = sum(ele_coords,2)/size( NodeIndices,2);

    % T0 = T0T1(:,1);
    % T1 = T0T1(:,2);
    % theta = VAT_fiber_ply_angle_1D(T0,T1,centerX',center,physical_length);

    c=cosd(theta);
    s=sind(theta);

    %%
    Q11b=Q11.*c.^4+2*(Q12+2.*Q66).*s.^2.*c.^2+Q22.*s.^4;
    Q12b=(Q11+Q22-4*Q66).*s.^2.*c.^2+Q12*(s.^4+c.^4);
    Q22b=Q11.*s.^4+2.*(Q12+2*Q66).*s.^2.*c.^2+Q22.*c.^4;
    Q16b=(Q11-Q12-2.*Q66).*s.*c.^3+(Q12-Q22+2*Q66).*s.^3.*c;
    Q26b=(Q11-Q12-2.*Q66).*c.*s.^3+(Q12-Q22+2*Q66).*c.^3.*s;
    Q66b=(Q11+Q22-2.*Q12-2.*Q66).*s.^2.*c.^2+Q66.*(s.^4+c.^4);
    % Q44b=Q44.*c.^2+Q55.*s.^2;
    % Q45b=(Q55-Q44).*c.*s;
    % Q55b=Q55.*c.^2+Q44.*s.^2;
            
    Q_mat=[Q11b,Q12b,Q16b,Q26b,Q66b,Q22b];