function Kplate = kstiff_assemble_3(FOM_ROM_flag,Stru,Laminate,Mat,FEM,T0T1,mat_files_folder)
    
    if strcmp(FOM_ROM_flag,'ROM') || strcmp(FOM_ROM_flag,'kROM_gFOM')
    
        center =Stru.center;
        physical_length = Stru.length; 
    
        % VAT fibers
        %NodeIndices=FEM.elementNodes(elem,:);%% Node NO. for one element
        
        elementNodes = FEM.elementNodes';
        elementNodes=reshape(elementNodes,[],1);
        ele_node_coord2 = FEM.nodeCoordinates_label(elementNodes,2);
        ele_node_coord2 = reshape(ele_node_coord2,size(FEM.elementNodes,2),[]);

        Element_center_X   = sum(ele_node_coord2)/size( FEM.elementNodes,2);

        K_mats = load([mat_files_folder filesep 'red_vat_elewise_mats.mat'],"RED_Kmn");
        % K_mats_red_reshp = structfun(@reshape_struct,K_mats.RED_Kmn,'uniformoutput',false);
        K_mats_red_reshp = K_mats.RED_Kmn;

        [V_mat,~]=lamination_parameters_elewise_array(Stru.thickness,Laminate.thickness,T0T1,Element_center_X,center,physical_length,Mat);
        % r_buck = size(K_mats.RED_Kmn.K_A0_red,1);
        V_mat_array =structfun(@(x) vat_reshape(x),V_mat,'uniformoutput',false);
        
        % size(V_mat_array.V0A)
        % size(K_mats_red_reshp.K_A0_red)
        Kplate = V_mat_array.V0A.*K_mats_red_reshp.K_A0_red  + V_mat_array.V1A.*K_mats_red_reshp.K_A1_red  + V_mat_array.V2A.*K_mats_red_reshp.K_A2_red  + ...
                    V_mat_array.V3A.*K_mats_red_reshp.K_A3_red  + V_mat_array.V4A.*K_mats_red_reshp.K_A4_red  +...
                 V_mat_array.V0B.*K_mats_red_reshp.K_B0_red  + V_mat_array.V1B.*K_mats_red_reshp.K_B1_red  + V_mat_array.V2B.*K_mats_red_reshp.K_B2_red  + ...
                    V_mat_array.V3B.*K_mats_red_reshp.K_B3_red  + V_mat_array.V4B.*K_mats_red_reshp.K_B4_red  +...
                 V_mat_array.V0D.*K_mats_red_reshp.K_D0_red  + V_mat_array.V1D.*K_mats_red_reshp.K_D1_red  + V_mat_array.V2D.*K_mats_red_reshp.K_D2_red  + ...
                    V_mat_array.V3D.*K_mats_red_reshp.K_D3_red + V_mat_array.V4D.*K_mats_red_reshp.K_D4_red +...
                 V_mat_array.V0A.*K_mats_red_reshp.K_Ash0_red + V_mat_array.V1A.*K_mats_red_reshp.K_Ash1_red + V_mat_array.V2A.*K_mats_red_reshp.K_Ash2_red;
% '-----'

        % size(Kplate)
        %         spy(Kplate)
                Kplate = reshape(full(Kplate),size(Kplate,1),size(Kplate,1),[]);
                % size(Kplate)
                Kplate = sum(Kplate,3);
                % size(Kplate)
                % Kplate_actv = Kplate;
        %         figure
        % spy(Kplate)

        % %ivec=reshape((1:r_buck)',1,[]);
        % ivec=repmat(1:r_buck,r_buck,FEM.elementNumber);
        % 
        % %jvec=repmat(elemdof_list,1,K_mat_size);
        % jvec=repmat((1:r_buck)',1,r_buck*FEM.elementNumber);
        % 
        % Sp_K_FOM.K_A0_3Dmat = sparse(ivec,jvec,V_mat_array.V0A.*K_mats_red_reshp.K_A0_red,r_buck,r_buck);
        % Sp_K_FOM.K_A1_3Dmat = sparse(ivec,jvec,V_mat_array.V1A.*K_mats_red_reshp.K_A1_red,r_buck,r_buck);
        % Sp_K_FOM.K_A2_3Dmat = sparse(ivec,jvec,V_mat_array.V2A.*K_mats_red_reshp.K_A2_red,r_buck,r_buck);
        % Sp_K_FOM.K_A3_3Dmat = sparse(ivec,jvec,V_mat_array.V3A.*K_mats_red_reshp.K_A3_red,r_buck,r_buck);
        % Sp_K_FOM.K_A4_3Dmat = sparse(ivec,jvec,V_mat_array.V4A.*K_mats_red_reshp.K_A4_red,r_buck,r_buck);
        % Sp_K_FOM.K_B0_3Dmat = sparse(ivec,jvec,V_mat_array.V0B.*K_mats_red_reshp.K_B0_red,r_buck,r_buck);
        % Sp_K_FOM.K_B1_3Dmat = sparse(ivec,jvec,V_mat_array.V1B.*K_mats_red_reshp.K_B1_red,r_buck,r_buck);
        % Sp_K_FOM.K_B2_3Dmat = sparse(ivec,jvec,V_mat_array.V2B.*K_mats_red_reshp.K_B2_red,r_buck,r_buck);
        % Sp_K_FOM.K_B3_3Dmat = sparse(ivec,jvec,V_mat_array.V3B.*K_mats_red_reshp.K_B3_red,r_buck,r_buck);
        % Sp_K_FOM.K_B4_3Dmat = sparse(ivec,jvec,V_mat_array.V4B.*K_mats_red_reshp.K_B4_red,r_buck,r_buck);
        % Sp_K_FOM.K_D0_3Dmat = sparse(ivec,jvec,V_mat_array.V0D.*K_mats_red_reshp.K_D0_red,r_buck,r_buck);
        % Sp_K_FOM.K_D1_3Dmat = sparse(ivec,jvec,V_mat_array.V1D.*K_mats_red_reshp.K_D1_red,r_buck,r_buck);
        % Sp_K_FOM.K_D2_3Dmat = sparse(ivec,jvec,V_mat_array.V2D.*K_mats_red_reshp.K_D2_red,r_buck,r_buck);
        % Sp_K_FOM.K_D3_3Dmat = sparse(ivec,jvec,V_mat_array.V3D.*K_mats_red_reshp.K_D3_red,r_buck,r_buck);
        % Sp_K_FOM.K_D4_3Dmat = sparse(ivec,jvec,V_mat_array.V4D.*K_mats_red_reshp.K_D4_red,r_buck,r_buck);
        % Sp_K_FOM.K_Ash0_3Dmat = sparse(ivec,jvec,V_mat_array.V2D.*K_mats_red_reshp.K_Ash0_red,r_buck,r_buck);
        % Sp_K_FOM.K_Ash1_3Dmat = sparse(ivec,jvec,V_mat_array.V3D.*K_mats_red_reshp.K_Ash1_red,r_buck,r_buck);
        % Sp_K_FOM.K_Ash2_3Dmat = sparse(ivec,jvec,V_mat_array.V4D.*K_mats_red_reshp.K_Ash2_red,r_buck,r_buck);
        % 
        % Kplate =Sp_K_FOM.K_A0_3Dmat + Sp_K_FOM.K_A1_3Dmat + Sp_K_FOM.K_A2_3Dmat   + ...
        %         Sp_K_FOM.K_A3_3Dmat + Sp_K_FOM.K_A4_3Dmat   +...
        %      Sp_K_FOM.K_B0_3Dmat + Sp_K_FOM.K_B1_3Dmat + Sp_K_FOM.K_B2_3Dmat   + ...
        %         Sp_K_FOM.K_B3_3Dmat + Sp_K_FOM.K_B4_3Dmat   +...
        %      Sp_K_FOM.K_D0_3Dmat + Sp_K_FOM.K_D1_3Dmat + Sp_K_FOM.K_D2_3Dmat   + ...
        %         Sp_K_FOM.K_D3_3Dmat  + Sp_K_FOM.K_D4_3Dmat  +...
        %      Sp_K_FOM.K_Ash0_3Dmat  + Sp_K_FOM.K_Ash1_3Dmat  + Sp_K_FOM.K_Ash2_3Dmat  ;

        % size(Kplate_fom)
        % Kplate_fom = reshape(Kplate_fom,size(Kplate_fom,1),size(Kplate_fom,1),[]);
        % size(Kplate_fom)
        % Kplate = sum(Kplate_fom,3);


        
        
        % row_indx = repmat([1:size(Kplate,1)]',1,size(Kplate,2)*FEM.numberElements);
        % col_indx = repmat(1:size(Kplate,2),size(Kplate,1),FEM.numberElements);
        % Kplate_actv = sparse(repmat([1:size(Kplate,1)]',1,size(Kplate,2)*FEM.numberElements),repmat(1:size(Kplate,2),size(Kplate,1),FEM.numberElements),Kplate);
        %Kplate = sum(Kplate,3);
        
        %Kplate_actv = Kplate;
    
    elseif strcmp(FOM_ROM_flag,'FOM')
    
        % The fiber path orientation is generalized for curvilinear fiber
        % determine the middle point in fiber path angle computation
        center =Stru.center;
        % length of the plate
        physical_length = Stru.length; % change along width direction
    
        [Kplate,Kelemp] = LinearStiffnessLaminatedPlate_VAT_v2_X_constant_angle(Mat, Stru,FEM,Laminate,T0T1,center(1),physical_length);
        
    end

end



function reshape_struct = reshape_struct(struct_var) %(struct_var,rsh_size)

    dim_var = size(struct_var);

    rsh_size=[dim_var(1),dim_var(2)*dim_var(3)];
    reshape_struct=reshape(struct_var,rsh_size);

end


function vat_struct = vat_reshape(struct_var)

    % u_r = load(['mat_files_latin_vat_stat280_buck91' filesep 'u_vat_9990.mat'],"r");
    % u_r = load(['mat_files_latin_vat' filesep 'case_buck117_stat155' filesep 'basis_vecs.mat'],"r_buck");
    % u_r = load(['mat_files_latin_vat' filesep 'case_n_train6400_buck32_stat47' filesep 'basis_vecs.mat'],"r_buck");
    u_r = load(['mat_files_latin_vat' filesep 'case_n_train400_buck24_stat42' filesep 'basis_vecs.mat'],"r_buck");

    %size(struct_var)
    B=sum(struct_var,1);
    %size(B)
    B=repmat(B',[1,u_r.r_buck]);
    %size(B)
    vat_struct=reshape(B',1,[]);
    %size(B)
    %vat_struct=sum(B,1);

end

