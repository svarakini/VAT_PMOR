function [KGplate_actv,FEM] = gstiff_assemble_rom_LS_dispfrmK_v2(FOM_ROM_flag,Stru,Laminate,Mat, FEM, T0T1,Kplate,mat_files_folder)

DisplacemntBoundaryCase=FEM.DisplacemntBoundaryCase;
center =Stru.center;
physical_length = Stru.length; 

if strcmp(FOM_ROM_flag,'ROM')
    
        %%%%%%%%%%%%%%%% Compute static displacement %%%%%%%%%%%%%%%%
    
        FEM.BCtype='SSSS-disp-Coburn-mid'; %% 'SSSS-NYY'; 'SSSS-NXXNYY'; 'SSSS-NXY'
        
        [ActiveDof,Constrained_Dof]=InPlaneBCPlate5Dof_v2(FEM);
        
        xx=FEM.nodeCoordinates(:,1);
        yy=FEM.nodeCoordinates(:,2);
        
        nodeNum=size(FEM.nodeCoordinates,1);
        SideNodesNumList_RHS = find(xx==Stru.length)';
        SideNodesNumList_LHS = find(xx==0)';
        
        % specify boundary conditions x=a, u = -1e-4;
        Dof_fixed_disp_LHS =  SideNodesNumList_LHS + 3*nodeNum ;
        Dof_fixed_disp_RHS =  SideNodesNumList_RHS + 3*nodeNum;
        
        ActiveDof_not_fixed  = setdiff(ActiveDof,[Dof_fixed_disp_LHS Dof_fixed_disp_RHS]);
        
        %%$ LOAD: Displacement Boundary Condition 
        %%% Options: 
        %%% 1. uniform 
        %%% 2. linear-Symm
        %%% 3. linear-Asymm
        %%% 4. linear-Symm-MidPoint-Low
        %%% 5. linear-Symm-MidPoint-High
        %%% 6. linear-Asymm-MidPoint-LHS:High-RHS:Low
        %%% 7. random
        
        % DisplacemntBoundaryCase='uniform'
        disp_mag=.1e4;
        u_limit_lhs=1; l_limit_lhs=0; % keep u_limit at 1
        u_limit_rhs=1; l_limit_rhs=0;
        
        nodenum_lhs=length(Dof_fixed_disp_LHS);
        nodenum_rhs=length(Dof_fixed_disp_RHS);
        
        [loaddistr_lhs,loaddistr_rhs]=displacementDistribution(DisplacemntBoundaryCase, ...
                                                                u_limit_lhs,l_limit_lhs,nodenum_lhs,...
                                                                u_limit_rhs,l_limit_rhs,nodenum_rhs);
        
        fixed_u =[ disp_mag*loaddistr_lhs' ; -disp_mag*loaddistr_rhs' ];

        %%
        K_mats_fom = load([mat_files_folder filesep 'red_vat_elewise_mats.mat'], 'FOM_K');

        elementNodes = FEM.elementNodes';
        elementNodes = reshape(elementNodes,[],1);
        ele_node_coord2 = FEM.nodeCoordinates_label(elementNodes,2);
        ele_node_coord2 = reshape(ele_node_coord2,size(FEM.elementNodes,2),[]);

        Element_center_X   = sum(ele_node_coord2)/size( FEM.elementNodes,2);

        [V_mat,Q_mat]=lamination_parameters_elewise_array(Stru.thickness,Laminate.thickness,T0T1,Element_center_X,center,physical_length,Mat);
        V_mat_array = structfun(@vat_reshape,V_mat,'uniformoutput',false);

        %K_mats_fom_reshp =structfun(@reshape_struct,K_mats_fom.FOM_K,'uniformoutput',false);
        K_mats_fom_reshp = K_mats_fom.FOM_K;
        K_mat_size = size(K_mats_fom_reshp.K_A0_3Dmat,1);

        ele_d = load([mat_files_folder filesep 'elemdof_list.mat']);
        elemdof_list = ele_d.elemdof_list;

        ivec=reshape(elemdof_list',1,[]);
        ivec=repmat(ivec,K_mat_size,1);

        jvec=repmat(elemdof_list,1,K_mat_size);
        jvec=reshape(jvec',K_mat_size,[]);

        size_Kfom = FEM.GDof;
        Sp_K_FOM.K_A0_3Dmat = sparse(ivec,jvec,V_mat_array.V0A.*K_mats_fom_reshp.K_A0_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_A1_3Dmat = sparse(ivec,jvec,V_mat_array.V1A.*K_mats_fom_reshp.K_A1_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_A2_3Dmat = sparse(ivec,jvec,V_mat_array.V2A.*K_mats_fom_reshp.K_A2_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_A3_3Dmat = sparse(ivec,jvec,V_mat_array.V3A.*K_mats_fom_reshp.K_A3_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_A4_3Dmat = sparse(ivec,jvec,V_mat_array.V4A.*K_mats_fom_reshp.K_A4_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B0_3Dmat = sparse(ivec,jvec,V_mat_array.V0B.*K_mats_fom_reshp.K_B0_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B1_3Dmat = sparse(ivec,jvec,V_mat_array.V1B.*K_mats_fom_reshp.K_B1_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B2_3Dmat = sparse(ivec,jvec,V_mat_array.V2B.*K_mats_fom_reshp.K_B2_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B3_3Dmat = sparse(ivec,jvec,V_mat_array.V3B.*K_mats_fom_reshp.K_B3_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B4_3Dmat = sparse(ivec,jvec,V_mat_array.V4B.*K_mats_fom_reshp.K_B4_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D0_3Dmat = sparse(ivec,jvec,V_mat_array.V0D.*K_mats_fom_reshp.K_D0_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D1_3Dmat = sparse(ivec,jvec,V_mat_array.V1D.*K_mats_fom_reshp.K_D1_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D2_3Dmat = sparse(ivec,jvec,V_mat_array.V2D.*K_mats_fom_reshp.K_D2_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D3_3Dmat = sparse(ivec,jvec,V_mat_array.V3D.*K_mats_fom_reshp.K_D3_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D4_3Dmat = sparse(ivec,jvec,V_mat_array.V4D.*K_mats_fom_reshp.K_D4_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_Ash0_3Dmat = sparse(ivec,jvec,V_mat_array.V2D.*K_mats_fom_reshp.K_Ash0_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_Ash1_3Dmat = sparse(ivec,jvec,V_mat_array.V3D.*K_mats_fom_reshp.K_Ash1_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_Ash2_3Dmat = sparse(ivec,jvec,V_mat_array.V4D.*K_mats_fom_reshp.K_Ash2_3Dmat,size_Kfom,size_Kfom);

        Kplate_fom = Sp_K_FOM.K_A0_3Dmat + Sp_K_FOM.K_A1_3Dmat + Sp_K_FOM.K_A2_3Dmat   + ...
                        Sp_K_FOM.K_A3_3Dmat + Sp_K_FOM.K_A4_3Dmat   +...
                     Sp_K_FOM.K_B0_3Dmat + Sp_K_FOM.K_B1_3Dmat + Sp_K_FOM.K_B2_3Dmat   + ...
                        Sp_K_FOM.K_B3_3Dmat + Sp_K_FOM.K_B4_3Dmat   +...
                     Sp_K_FOM.K_D0_3Dmat + Sp_K_FOM.K_D1_3Dmat + Sp_K_FOM.K_D2_3Dmat   + ...
                        Sp_K_FOM.K_D3_3Dmat  + Sp_K_FOM.K_D4_3Dmat  +...
                     Sp_K_FOM.K_Ash0_3Dmat  + Sp_K_FOM.K_Ash1_3Dmat  + Sp_K_FOM.K_Ash2_3Dmat  ;

        K21 = Kplate_fom( ActiveDof_not_fixed ,  [Dof_fixed_disp_LHS Dof_fixed_disp_RHS]);
        K22 = Kplate_fom( ActiveDof_not_fixed ,  ActiveDof_not_fixed);
        
        % final_unknowns = -K22\(K21*fixed_u);
               
        % discrete analytical stiffness matrix
        % [ K11 K12 ]{X1}  {F1}
        % [         ]    =
        % [ K21 K22 ]{X2}  {0}
        %
        % K21 X1 + K22 X2 = 0; BASED ON X1 -> X2
        
        %%
        % FEM.displacement_fom=zeros(1,FEM.GDof);
        % FEM.displacement_fom( [Dof_fixed_disp_LHS Dof_fixed_disp_RHS] ) = fixed_u';
        % FEM.displacement_fom( ActiveDof_not_fixed) =  final_unknowns';

        Vbasismat=load([mat_files_folder filesep 'basis_vecs.mat']);
        % Vde = Vbasismat.Vde;
        r_buck = Vbasismat.r_buck;
        r_stat = Vbasismat.r_stat;

        % % 
        % % Kplate_fom=zeros(FEM.GDof);
        % % for elem=1:FEM.numberElements
        % %     NodeIndices=FEM.elementNodes(elem,:);%% Node NO. for one element
        % %     Element_center_X   = sum(FEM.nodeCoordinates_label( NodeIndices,2))/length( NodeIndices);
        % %     V_mat=lamination_parameters_elewise(Stru.thickness,Laminate.thickness,T0T1,Element_center_X,center,physical_length);
        % %     dof_indx=elemdof_list(elem,:);
        % %     Kplate_fom(dof_indx,dof_indx) = Kplate_fom(dof_indx,dof_indx) + ...
        % %          sum(V_mat.V0A)*K_mats_fom.FOM_K.K_A0_3Dmat(:,:,elem) + sum(V_mat.V1A)*K_mats_fom.FOM_K.K_A1_3Dmat(:,:,elem) + sum(V_mat.V2A)*K_mats_fom.FOM_K.K_A2_3Dmat(:,:,elem) + ...
        % %             sum(V_mat.V3A)*K_mats_fom.FOM_K.K_A3_3Dmat(:,:,elem) + sum(V_mat.V4A)*K_mats_fom.FOM_K.K_A4_3Dmat(:,:,elem) +...
        % %          sum(V_mat.V0B)*K_mats_fom.FOM_K.K_B0_3Dmat(:,:,elem) + sum(V_mat.V1B)*K_mats_fom.FOM_K.K_B1_3Dmat(:,:,elem) + sum(V_mat.V2B)*K_mats_fom.FOM_K.K_B2_3Dmat(:,:,elem) + ...
        % %             sum(V_mat.V3B)*K_mats_fom.FOM_K.K_B3_3Dmat(:,:,elem) + sum(V_mat.V4B)*K_mats_fom.FOM_K.K_B4_3Dmat(:,:,elem) +...
        % %          sum(V_mat.V0D)*K_mats_fom.FOM_K.K_D0_3Dmat(:,:,elem) + sum(V_mat.V1D)*K_mats_fom.FOM_K.K_D1_3Dmat(:,:,elem) + sum(V_mat.V2D)*K_mats_fom.FOM_K.K_D2_3Dmat(:,:,elem) + ...
        % %             sum(V_mat.V3D)*K_mats_fom.FOM_K.K_D3_3Dmat(:,:,elem)+ sum(V_mat.V4D)*K_mats_fom.FOM_K.K_D4_3Dmat(:,:,elem)+...
        % %          sum(V_mat.V0A)*K_mats_fom.FOM_K.K_Ash0_3Dmat(:,:,elem)+ sum(V_mat.V1A)*K_mats_fom.FOM_K.K_Ash1_3Dmat(:,:,elem)+ sum(V_mat.V2A)*K_mats_fom.FOM_K.K_Ash2_3Dmat(:,:,elem);
        % % end
        % fixed_dof=[Dof_fixed_disp_LHS Dof_fixed_disp_RHS];
        % % % % % % % lhs = Vbasismat.Vde;
        % % % % % % % 
        % % % % % % % % [~,~,pttss]=qr(lhs','vector');
        % % % % % % % ppttss=load([mat_files_folder filesep 'pts_test']);
        % % % % % % % pts=ppttss.pts(1:r_stat);
        % % % % % % % 
        % % % % % % % rhs = final_unknowns; % FEM.displacement_fom(ActiveDof_not_fixed);
        % % % % % % % 
        % % % % % % % %%
        % % % % % % % rhs = K21*fixed_u;
        % % % % % % % % % lhs = K21(pts,:)*fixed_u;
        % % % % % % % lhs = -K22(:,pts)*Vbasismat.Vde(pts,:); %final_unknowns;


        rhs = K21*fixed_u;
        % % lhs = K21(pts,:)*fixed_u;
        lhs = -K22*Vbasismat.Vde;
        %%
        % %         lhs_rhs=load([mat_files_folder filesep 'Krylov_vecs_04182023.mat']);

        % %         pts=Vbasismat.pts;
        % %         elementNodes = FEM.elementNodes';
        % %         elementNodes=reshape(elementNodes,[],1);
        % %         ele_node_coord2 = FEM.nodeCoordinates_label(elementNodes,2);
        % %         ele_node_coord2 = reshape(ele_node_coord2,size(FEM.elementNodes,2),[]);
        % %         
        % %         Element_center_X   = sum(ele_node_coord2)/size( FEM.elementNodes,2);
        % % 
        % %         V_mat = lamination_parameters_elewise_array(Stru.thickness,Laminate.thickness,T0T1,Element_center_X,center,physical_length);
        % %         V_laminates_lhs = structfun(@vat_reshape_lhs,V_mat,'uniformoutput',false);
        % % 
        % %         lhs_vecs = V_laminates_lhs.V0A.*lhs_rhs.lhs_K_A0_3Dmat  + V_laminates_lhs.V1A.*lhs_rhs.lhs_K_A1_3Dmat  + V_laminates_lhs.V2A.*lhs_rhs.lhs_K_A2_3Dmat  + ...
        % %                     V_laminates_lhs.V3A.*lhs_rhs.lhs_K_A3_3Dmat  + V_laminates_lhs.V4A.*lhs_rhs.lhs_K_A4_3Dmat  +...
        % %                  V_laminates_lhs.V0B.*lhs_rhs.lhs_K_B0_3Dmat  + V_laminates_lhs.V1B.*lhs_rhs.lhs_K_B1_3Dmat  + V_laminates_lhs.V2B.*lhs_rhs.lhs_K_B2_3Dmat  + ...
        % %                     V_laminates_lhs.V3B.*lhs_rhs.lhs_K_B3_3Dmat  + V_laminates_lhs.V4B.*lhs_rhs.lhs_K_B4_3Dmat  +...
        % %                  V_laminates_lhs.V0D.*lhs_rhs.lhs_K_D0_3Dmat  + V_laminates_lhs.V1D.*lhs_rhs.lhs_K_D1_3Dmat  + V_laminates_lhs.V2D.*lhs_rhs.lhs_K_D2_3Dmat  + ...
        % %                     V_laminates_lhs.V3D.*lhs_rhs.lhs_K_D3_3Dmat + V_laminates_lhs.V4D.*lhs_rhs.lhs_K_D4_3Dmat +...
        % %                  V_laminates_lhs.V0A.*lhs_rhs.lhs_K_Ash0_3Dmat + V_laminates_lhs.V1A.*lhs_rhs.lhs_K_Ash1_3Dmat + V_laminates_lhs.V2A.*lhs_rhs.lhs_K_Ash2_3Dmat ;
        % %         lhs = -sum(lhs_vecs,3);
        % % 
        % %         V_laminates_rhs =structfun(@vat_reshape_rhs,V_mat,'uniformoutput',false);
        % %         rhs_vecs = lhs_rhs.rhs; 
        % %         rhs_vec_assebled = V_laminates_rhs.V0A.*rhs_vecs(:,:,1)  + V_laminates_rhs.V1A.*rhs_vecs(:,:,2)  + V_laminates_rhs.V2A.*rhs_vecs(:,:,3)  + ...
        % %                     V_laminates_rhs.V3A.*rhs_vecs(:,:,4)  + V_laminates_rhs.V4A.*rhs_vecs(:,:,5)  +...
        % %                  V_laminates_rhs.V0B.*rhs_vecs(:,:,6)  + V_laminates_rhs.V1B.*rhs_vecs(:,:,7)  + V_laminates_rhs.V2B.*rhs_vecs(:,:,8)  + ...
        % %                     V_laminates_rhs.V3B.* rhs_vecs(:,:,9)  + V_laminates_rhs.V4B.* rhs_vecs(:,:,10)  +...
        % %                  V_laminates_rhs.V0D.* rhs_vecs(:,:,11)  + V_laminates_rhs.V1D.* rhs_vecs(:,:,12)  + V_laminates_rhs.V2D.* rhs_vecs(:,:,13)  + ...
        % %                     V_laminates_rhs.V3D.* rhs_vecs(:,:,14) + V_laminates_rhs.V4D.* rhs_vecs(:,:,15) +...
        % %                  V_laminates_rhs.V0A.* rhs_vecs(:,:,16) + V_laminates_rhs.V1A.* rhs_vecs(:,:,17) + V_laminates_rhs.V2A.* rhs_vecs(:,:,18) ;
        % %         
        % %         rhs = sum(rhs_vec_assebled,2);
        
        %%
        % coefs = lhs(pts,:)\rhs(pts);

        coefs = lhs\rhs;
        %save([pwd filesep 'mat_files' filesep 'coef_vals.mat'],"coefs")

        % vde_mat=load([pwd filesep 'mat_files' filesep 'Vde_mat.mat']);
        % vde_mat=load([pwd filesep 'mat_files' filesep 'Vde_mat_653.mat']);

        % size(coefs)
        % size(vde_mat.Vde)
        % FEM.displacement =zeros(1,FEM.GDof);
        % FEM.displacement(ActiveDof_not_fixed) = Vbasismat.Vde*coefs; 
        % FEM.displacement([Dof_fixed_disp_LHS Dof_fixed_disp_RHS]) = fixed_u; 
        %size_u_mixed= size(u_mixed,2)
        % r_stat = size(coefs,1);

        % cc = reshape(coefs',[1,1,r_stat]);       

      % size(coefs)
        % size(vde_mat.Vde)
        % FEM.displacement =zeros(1,FEM.GDof);
        % FEM.displacement(ActiveDof_not_fixed) = Vbasismat.Vde*coefs; 
        % FEM.displacement([Dof_fixed_disp_LHS Dof_fixed_disp_RHS]) = fixed_u; 
        %size_u_mixed= size(u_mixed,2)
        

        cc = reshape(coefs',[1,1,r_stat]);       

        % G_mat = load([mat_files_folder filesep 'G_stiff_st_u_krylov_latin_r280']);
        stress_mat = load([mat_files_folder filesep 'stress_mats.mat']);
        %G_mat = load([pwd filesep 'mat_files' filesep 'G_stiff_vat_u_krylov_05222024']);
        % G_mat = load([pwd filesep 'mat_files' filesep 'G_stiff_vat_u_krylov_r653']);
        % stress_mat = load([pwd filesep 'mat_files' filesep 'stress_mats.mat']);
        % stress_mat = load([pwd filesep 'mat_files' filesep 'stress_mats_r653.mat']);
        
        Q_mat_n_compnts = 6;
        
        zcoordinate = -Stru.thickness/2 +Laminate.thickness(1:size(T0T1,1))/2 +(0:size(T0T1,1)-1).*Laminate.thickness;

        QQ_ele = permute(reshape(Q_mat',[FEM.numberElements,Q_mat_n_compnts,Laminate.layer,1]),[1,4,3,2]);
        %QQ_curv_ele = reshape(sum(Q_mat),[FEM.numberElements,1,1,6]);
        z_rshp = reshape(zcoordinate,[1,1,Laminate.layer]);

        FEM.stress = sum(QQ_ele.*repmat(sum(cc.*stress_mat.stress_normal_mat(:,:,1:r_stat,:),3),[1,1,Laminate.layer,1]) ...
                        + QQ_ele.*z_rshp.*repmat(sum(cc.*stress_mat.stress_curv_mat(:,:,1:r_stat,:),3),[1,1,Laminate.layer,1]) ...
                        + QQ_ele.*repmat(stress_mat.stress_normal_mat(:,:,end,:),[1,1,Laminate.layer,1]) ...
                        + QQ_ele.*z_rshp.*repmat(stress_mat.stress_curv_mat(:,:,end,:),[1,1,Laminate.layer,1]),4);

        layer=1:8;
        zbot = -Stru.thickness/2  + (layer-1).*Laminate.thickness(layer);
        ztop = -Stru.thickness/2  + layer.*Laminate.thickness(layer);
      
        G_mat = load([mat_files_folder filesep 'G_red_mats.mat']);
        G_red = G_mat.G_red;
        GG = reshape(G_red,size(G_red,1),[]);

        hhh=reshape(pagetranspose(repmat(reshape(FEM.stress,FEM.numberElements*3,1,[]),1,r_buck)),1,[],size(FEM.stress,3));
        KGplate_actv = sum(reshape(hhh.*reshape(ztop-zbot,1,1,[]).*repmat(GG,1,1,8),r_buck,r_buck,[]),3);



  
elseif strcmp(FOM_ROM_flag,'FOM')
    
        %%%%%%%%%% STATIC ANALYSIS AT INITIAL DISPLACEMENT %%%%%%%%%
        
        %%================== Boundary Conditions ===================
        FEM.BCtype='SSSS-disp-Coburn-mid'; %% 'SSSS-NYY'; 'SSSS-NXXNYY'; 'SSSS-NXY'
        
        [ActiveDof,Constrained_Dof]=InPlaneBCPlate5Dof_v2(FEM);
        
        xx=FEM.nodeCoordinates(:,1);
        yy=FEM.nodeCoordinates(:,2);
        
        nodeNum=size(FEM.nodeCoordinates,1);
        SideNodesNumList_RHS = find(xx==Stru.length)';
        SideNodesNumList_LHS = find(xx==0)';
        
        % specify boundary conditions x=a, u = -1e-4;
        Dof_fixed_disp_LHS =  SideNodesNumList_LHS + 3*nodeNum ;
        Dof_fixed_disp_RHS =  SideNodesNumList_RHS + 3*nodeNum;
        
        ActiveDof_not_fixed  = setdiff(ActiveDof,[Dof_fixed_disp_LHS Dof_fixed_disp_RHS]);
        
        %%============= LOAD: Displacement Boundary Condition ==============
        %%% Options: 
        %%% 1. uniform 
        %%% 2. linear-Symm
        %%% 3. linear-Asymm
        %%% 4. linear-Symm-MidPoint-Low
        %%% 5. linear-Symm-MidPoint-High
        %%% 6. linear-Asymm-MidPoint-LHS:High-RHS:Low
        %%% 7. random
        
        % DisplacemntBoundaryCase='uniform'
        disp_mag=.1e4;
        u_limit_lhs=1; l_limit_lhs=0; % keep u_limit at 1
        u_limit_rhs=1; l_limit_rhs=0;
        
        nodenum_lhs=length(Dof_fixed_disp_LHS);
        nodenum_rhs=length(Dof_fixed_disp_RHS);
        
        [loaddistr_lhs,loaddistr_rhs]=displacementDistribution(DisplacemntBoundaryCase, ...
                                                                u_limit_lhs,l_limit_lhs,nodenum_lhs,...
                                                                u_limit_rhs,l_limit_rhs,nodenum_rhs);
        
        fixed_u =[ disp_mag*loaddistr_lhs' ; -disp_mag*loaddistr_rhs' ];
        
        K21 = Kplate( ActiveDof_not_fixed ,  [Dof_fixed_disp_LHS Dof_fixed_disp_RHS]);
        K22 = Kplate( ActiveDof_not_fixed ,  ActiveDof_not_fixed);
        
        final_unknowns = -K22\(K21*fixed_u);
        
        % discrete analytical stiffness matrix
        % [ K11 K12 ]{X1}  {F1}
        % [         ]    =
        % [ K21 K22 ]{X2}  {0}
        %
        % K21 X1 + K22 X2 = 0; BASED ON X1 -> X2
        
        FEM.displacement=zeros(1,FEM.GDof);
        
        FEM.displacement( [Dof_fixed_disp_LHS Dof_fixed_disp_RHS] ) = fixed_u';
        FEM.displacement( ActiveDof_not_fixed) =  final_unknowns';

        center =Stru.center;
        [FEM.stress,FEM.strain]=StressRecoveryPlate_VAT_center_average_v2_constant_angle(FEM,Laminate,Mat,Stru,T0T1,center,Stru.length);
    
        % GEOMETRIIC STIFFNESS
        KGplate = GeometryStiffnessPlate_stress_recovery_noshear(FEM,Stru,Laminate);
        KGplate_actv = KGplate(FEM.ActiveDof_ssss,FEM.ActiveDof_ssss);


elseif strcmp(FOM_ROM_flag,'kROM_gFOM')
    
        %%%%%%%%%%%%%%%% Compute static displacement %%%%%%%%%%%%%%%%
    
        FEM.BCtype='SSSS-disp-Coburn-mid'; %% 'SSSS-NYY'; 'SSSS-NXXNYY'; 'SSSS-NXY'
        
        [ActiveDof,~]=InPlaneBCPlate5Dof_v2(FEM);
        
        xx=FEM.nodeCoordinates(:,1);
        %yy=FEM.nodeCoordinates(:,2);
        
        nodeNum=size(FEM.nodeCoordinates,1);
        SideNodesNumList_RHS = find(xx==Stru.length)';
        SideNodesNumList_LHS = find(xx==0)';
        
        % specify boundary conditions x=a, u = -1e-4;
        Dof_fixed_disp_LHS =  SideNodesNumList_LHS + 3*nodeNum ;
        Dof_fixed_disp_RHS =  SideNodesNumList_RHS + 3*nodeNum;
        
        ActiveDof_not_fixed  = setdiff(ActiveDof,[Dof_fixed_disp_LHS Dof_fixed_disp_RHS]);
        
        %%$ LOAD: Displacement Boundary Condition 
        %%% Options: 
        %%% 1. uniform 
        %%% 2. linear-Symm
        %%% 3. linear-Asymm
        %%% 4. linear-Symm-MidPoint-Low
        %%% 5. linear-Symm-MidPoint-High
        %%% 6. linear-Asymm-MidPoint-LHS:High-RHS:Low
        %%% 7. random
        
        % DisplacemntBoundaryCase='uniform'
        disp_mag=1; %.1e4;
        u_limit_lhs=1; l_limit_lhs=0; % keep u_limit at 1
        u_limit_rhs=1; l_limit_rhs=0;
        
        nodenum_lhs=length(Dof_fixed_disp_LHS);
        nodenum_rhs=length(Dof_fixed_disp_RHS);
        
        [loaddistr_lhs,loaddistr_rhs]=displacementDistribution(DisplacemntBoundaryCase, ...
                                                                u_limit_lhs,l_limit_lhs,nodenum_lhs,...
                                                                u_limit_rhs,l_limit_rhs,nodenum_rhs);
        
        fixed_u =[ disp_mag*loaddistr_lhs' ; -disp_mag*loaddistr_rhs' ];
    
        % V_mat=lamination_parameters(Stru.thickness,Laminate.thickness,T0T1);
        K_mats_fom=load([mat_files_folder filesep 'red_vat_elewise_mats.mat'], 'FOM_K');

        elementNodes = FEM.elementNodes';
        elementNodes=reshape(elementNodes,[],1);
        ele_node_coord2 = FEM.nodeCoordinates_label(elementNodes,2);
        ele_node_coord2 = reshape(ele_node_coord2,size(FEM.elementNodes,2),[]);

        Element_center_X   = sum(ele_node_coord2)/size( FEM.elementNodes,2);

        V_mat=lamination_parameters_elewise_array(Stru.thickness,Laminate.thickness,T0T1,Element_center_X,center,physical_length,Mat);
        V_mat_array =structfun(@vat_reshape,V_mat,'uniformoutput',false);

        %K_mats_fom_reshp =structfun(@reshape_struct,K_mats_fom.FOM_K,'uniformoutput',false);
        K_mats_fom_reshp = K_mats_fom.FOM_K;
        K_mat_size = size(K_mats_fom_reshp.K_A0_3Dmat,1);


        ele_d = load([mat_files_folder filesep 'elemdof_list.mat']);
        elemdof_list = ele_d.elemdof_list;

        ivec=reshape(elemdof_list',1,[]);
        ivec=repmat(ivec,K_mat_size,1);

        jvec=repmat(elemdof_list,1,K_mat_size);
        jvec=reshape(jvec',K_mat_size,[]);

        
        size_Kfom = FEM.GDof;
        Sp_K_FOM.K_A0_3Dmat = sparse(ivec,jvec,V_mat_array.V0A.*K_mats_fom_reshp.K_A0_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_A1_3Dmat = sparse(ivec,jvec,V_mat_array.V1A.*K_mats_fom_reshp.K_A1_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_A2_3Dmat = sparse(ivec,jvec,V_mat_array.V2A.*K_mats_fom_reshp.K_A2_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_A3_3Dmat = sparse(ivec,jvec,V_mat_array.V3A.*K_mats_fom_reshp.K_A3_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_A4_3Dmat = sparse(ivec,jvec,V_mat_array.V4A.*K_mats_fom_reshp.K_A4_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B0_3Dmat = sparse(ivec,jvec,V_mat_array.V0B.*K_mats_fom_reshp.K_B0_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B1_3Dmat = sparse(ivec,jvec,V_mat_array.V1B.*K_mats_fom_reshp.K_B1_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B2_3Dmat = sparse(ivec,jvec,V_mat_array.V2B.*K_mats_fom_reshp.K_B2_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B3_3Dmat = sparse(ivec,jvec,V_mat_array.V3B.*K_mats_fom_reshp.K_B3_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_B4_3Dmat = sparse(ivec,jvec,V_mat_array.V4B.*K_mats_fom_reshp.K_B4_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D0_3Dmat = sparse(ivec,jvec,V_mat_array.V0D.*K_mats_fom_reshp.K_D0_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D1_3Dmat = sparse(ivec,jvec,V_mat_array.V1D.*K_mats_fom_reshp.K_D1_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D2_3Dmat = sparse(ivec,jvec,V_mat_array.V2D.*K_mats_fom_reshp.K_D2_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D3_3Dmat = sparse(ivec,jvec,V_mat_array.V3D.*K_mats_fom_reshp.K_D3_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_D4_3Dmat = sparse(ivec,jvec,V_mat_array.V4D.*K_mats_fom_reshp.K_D4_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_Ash0_3Dmat = sparse(ivec,jvec,V_mat_array.V2D.*K_mats_fom_reshp.K_Ash0_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_Ash1_3Dmat = sparse(ivec,jvec,V_mat_array.V3D.*K_mats_fom_reshp.K_Ash1_3Dmat,size_Kfom,size_Kfom);
        Sp_K_FOM.K_Ash2_3Dmat = sparse(ivec,jvec,V_mat_array.V4D.*K_mats_fom_reshp.K_Ash2_3Dmat,size_Kfom,size_Kfom);


        Kplate_fom =Sp_K_FOM.K_A0_3Dmat + Sp_K_FOM.K_A1_3Dmat + Sp_K_FOM.K_A2_3Dmat   + ...
                        Sp_K_FOM.K_A3_3Dmat + Sp_K_FOM.K_A4_3Dmat   +...
                     Sp_K_FOM.K_B0_3Dmat + Sp_K_FOM.K_B1_3Dmat + Sp_K_FOM.K_B2_3Dmat   + ...
                        Sp_K_FOM.K_B3_3Dmat + Sp_K_FOM.K_B4_3Dmat   +...
                     Sp_K_FOM.K_D0_3Dmat + Sp_K_FOM.K_D1_3Dmat + Sp_K_FOM.K_D2_3Dmat   + ...
                        Sp_K_FOM.K_D3_3Dmat  + Sp_K_FOM.K_D4_3Dmat  +...
                     Sp_K_FOM.K_Ash0_3Dmat  + Sp_K_FOM.K_Ash1_3Dmat  + Sp_K_FOM.K_Ash2_3Dmat  ;

        % Kplate_fom = reshape(Kplate_fom,size(Kplate_fom,1),size(Kplate_fom,1),[]);
        % Kplate_fom = sum(Kplate_fom,3);
        
        
        K21 = Kplate_fom( ActiveDof_not_fixed ,  [Dof_fixed_disp_LHS Dof_fixed_disp_RHS]);
        K22 = Kplate_fom( ActiveDof_not_fixed ,  ActiveDof_not_fixed);
        
        final_unknowns = -K22\(K21*fixed_u);
               
        % discrete analytical stiffness matrix
        % [ K11 K12 ]{X1}  {F1}
        % [         ]    =
        % [ K21 K22 ]{X2}  {0}
        %
        % K21 X1 + K22 X2 = 0; BASED ON X1 -> X2
        
        FEM.displacement=zeros(1,FEM.GDof);
        
        FEM.displacement( [Dof_fixed_disp_LHS Dof_fixed_disp_RHS] ) = fixed_u';
        FEM.displacement( ActiveDof_not_fixed) =  final_unknowns';
    
    
        [FEM.stress,FEM.strain]=StressRecoveryPlate_VAT_center_average_v2_constant_angle(FEM,Laminate,Mat,Stru,T0T1,center,Stru.length);
    
        % GEOMETRIIC STIFFNESS
        KGplate = GeometryStiffnessPlate_stress_recovery(FEM,Stru,Laminate);

        basis_mat=load([mat_files_folder filesep 'basis_vecs.mat'] );
        % basis_mat=load([mat_files_folder filesep 'u_vat_9990.mat'] );
        u_mixed=basis_mat.u_mixed;
        
        KGplate_actv = u_mixed'*KGplate(FEM.ActiveDof_ssss,FEM.ActiveDof_ssss)*u_mixed;       

end

end



function reshape_struct = reshape_struct(struct_var) %(struct_var,rsh_size)

    dim_var = size(struct_var);
    rsh_size=[dim_var(1),dim_var(2)*dim_var(3)];
    reshape_struct=reshape(struct_var,rsh_size);

end


function vat_struct = vat_reshape(struct_var) %(struct_var,rsh_size)

    %global FEM

    B=sum(struct_var,1);
    B=repmat(B',[1,40]);
    vat_struct=reshape(B',1,[]);

end

function vat_struct = vat_reshape_lhs(struct_var) %(struct_var,rsh_size)
    B=sum(struct_var,1);
    vat_struct=reshape(B,1,1,size(B,2));
end



function vat_struct = vat_reshape_rhs(struct_var) %(struct_var,rsh_size)
    B=sum(struct_var,1);
    vat_struct=B;
end
