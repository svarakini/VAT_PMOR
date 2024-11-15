% pic_folder = "C:\Users\svarakini\OneDrive - Virginia Tech\Desktop\1_Thesis_PlateMOR_2_07312023\1_new_from_arc_updated_May2024\K_affine_only\b1_Plate_PSO_elewiseROM\mode_plots";
% pic_folder="C:\Users\svarakini\Desktop\G_stiff_05222024\b1_g_rom_LS\b1_g_rom_LS\subroutines\platebuckling_unstiff_FEA_MOR\mode_plots\stress_plots\basis_vector\"
% plot_folder = ['..\..\mode_plots\stress_plots\test_2\basis_vector\'];

save_fig = 1;
plot_folder = ['..\..\mode_plots\trial1\'];


%%
% close all
ActiveDof=FEM.ActiveDof_ssss;

% file_name_mode = strcat("modeshape_m27p5_", FOM_ROM_flag); 
% file_name_modeerror = strcat("modeshape_error_m27p5_", FOM_ROM_flag); 
% file_name_mode = strcat("modeshape_m7p5_50_", FOM_ROM_flag); %fullfile("modeshape_m7p5_50_", FOM_ROM_flag);
% file_name_modeerror = strcat("modeshape_error_m7p5_50", FOM_ROM_flag); %fullfile("modeshape_error_m7p5_50", FOM_ROM_flag);
% file_name_mode = strcat("modeshape_QI_", FOM_ROM_flag); 
% file_name_modeerror = strcat("modeshape_error_QI_", FOM_ROM_flag); 


static_displacement=FEM.displacement';
% basis_num=4;
% static_displacement=Vde(:,basis_num); %FEM.displacement';
% FOM_ROM_flag="Vbasis"
% filename_mode = char(strcat("modeshape_", FOM_ROM_flag,num2str(basis_num))); 


% load(["mat_files\u_vat_test_9990.mat"])
% load(["mat_files\Vde_mat.mat"])
% if ~strcmp(FOM_ROM_flag,'FOM')
%     static_displacement=Vde*static_displacement';
% end

X=zeros(size(FEM.elementNodes,1),size(FEM.nodesCord,2));
Y=zeros(size(FEM.elementNodes,1),size(FEM.nodesCord,2));
Z=zeros(size(FEM.elementNodes,1),size(FEM.nodesCord,2));

coordinates=FEM.nodesCord;
nodes(1:size(FEM.elementNodes,1),1)=1:size(FEM.elementNodes,1);
nodes(1:FEM.elementNumber,2:size(FEM.elementNodes,2)+1)=FEM.elementNodes;
%         %---Plot Mesh----
% PlotMesh(coordinates,nodes)
% view(2)
%

deformUZ=zeros(FEM.NodeNumber,1);
deformBx=zeros(FEM.NodeNumber,1);
deformBy=zeros(FEM.NodeNumber,1);
deformUX=zeros(FEM.NodeNumber,1);
deformUY=zeros(FEM.NodeNumber,1);


bendingmode1=find(ActiveDof<=FEM.NodeNumber);
% bendingmode=find(ActiveDof>FEM.NodeNumber & ActiveDof<=2*FEM.NodeNumber );

ActiveBendDOF=ActiveDof(bendingmode1); %-FEM.NodeNumber;

bendModeNo=length(bendingmode1);

%---------------Mode Shape Plot-------------------

if bendModeNo<10
    mode=bendModeNo;
else
    mode=plotmodenNo;
end

% for modeNo=  mode
    
    deformUZ(ActiveBendDOF)=static_displacement(1:bendModeNo,:);
    dx=min(Xcoord):(max(Xcoord)-min(Xcoord))/100:max(Xcoord);
    dy=min(Ycoord):(max(Ycoord)-min(Ycoord))/100:max(Ycoord);
    [x3,y3]=meshgrid(dx,dy,0);
    z3=griddata(Xcoord,Ycoord,deformUZ,x3,y3,'v4');
    % ZQ=griddata(Xcoord,Ycoord,Zcoord,deformUZ,x3,y3,z3,'linear');
    hf=figure;
    axes1 = axes('Parent',hf,'YTick',[0 0.8 0.801],'XTick',[0 0.6 0.601],...
        'DataAspectRatio',[1 1 1]);
    hold(axes1,'all');
    
    [max_value,max_id]=max(abs(z3(:)));
    set(gcf,'color','w')
    
    % scalefactor = 1/(max(abs(z3(:)))/max([max(abs(y3(:))) max(abs(x3(:)))]))/Stru.width;
    % 
    % if max(abs(z3(:))) == max((z3(:)))
    %     sign_max = 1;
    % elseif max(abs(z3(:))) == -min((z3(:)))
    %     sign_max = -1;
    % end
    
    fig_mode = surf(x3,y3,z3,'FaceColor','interp',...
        'LineStyle','none',...
        'EdgeColor','none',...
        'FaceLighting','phong');hold on;
    % clim([-7.145577607161163e-12 5.803395986890626e-16])
    % clim([-6.426388468125317e-11 1.869523785304594e-10])
    % clim([-4.656606151909073e-11 1.651731625297348e-12])
    % clim( [-7.178057719384481e-12 3.405127415371306e-13])
    hh=colorbar('FontSize',16)
    % colorbar('FontSize',16,'Limits',[-3e-11,0.25e-11])
    colormap(jet(20));

    if save_fig == 1
    saveas(fig_mode,[plot_folder filename_mode '_z'],file_typ)
    end


%%
% bendingmode=find(ActiveDof<=FEM.NodeNumber);
bendingmode2=find(ActiveDof>FEM.NodeNumber & ActiveDof<=2*FEM.NodeNumber );

ActiveBendDOF=ActiveDof(bendingmode2)-FEM.NodeNumber;

bendModeNo=length(bendingmode2);

%---------------Mode Shape Plot-------------------

if bendModeNo<10
    mode=bendModeNo;
else
    mode=plotmodenNo;
end

% for modeNo=  mode
    
    deformBx(ActiveBendDOF)=static_displacement(1:bendModeNo,:);
    dx=min(Xcoord):(max(Xcoord)-min(Xcoord))/100:max(Xcoord);
    dy=min(Ycoord):(max(Ycoord)-min(Ycoord))/100:max(Ycoord);
    [x3,y3]=meshgrid(dx,dy,0);
    z3=griddata(Xcoord,Ycoord,deformBx,x3,y3,'v4');
    % ZQ=griddata(Xcoord,Ycoord,Zcoord,deformUZ,x3,y3,z3,'linear');
    hf=figure;
    axes1 = axes('Parent',hf,'YTick',[0 0.8 0.801],'XTick',[0 0.6 0.601],...
        'DataAspectRatio',[1 1 1]);
    hold(axes1,'all');
    
    [max_value,max_id]=max(abs(z3(:)));
    set(gcf,'color','w')
    
    % scalefactor = 1/(max(abs(z3(:)))/max([max(abs(y3(:))) max(abs(x3(:)))]))/Stru.width;
    % 
    % if max(abs(z3(:))) == max((z3(:)))
    %     sign_max = 1;
    % elseif max(abs(z3(:))) == -min((z3(:)))
    %     sign_max = -1;
    % end
    
    fig_mode = surf(x3,y3,z3,'FaceColor','interp',...
        'LineStyle','none',...
        'EdgeColor','none',...
        'FaceLighting','phong');hold on;
    % clim([-7.145577607161163e-12 5.803395986890626e-16])
    % clim([-6.426388468125317e-11 1.869523785304594e-10])
    % clim([-4.656606151909073e-11 1.651731625297348e-12])
    % clim( [-7.178057719384481e-12 3.405127415371306e-13])
    hh=colorbar('FontSize',16)
    % colorbar('FontSize',16,'Limits',[-3e-11,0.25e-11])
    colormap(jet(20));

    if save_fig == 1
    saveas(fig_mode,[plot_folder filename_mode '_Rx'],file_typ)
    end







%%
% bendingmode=find(ActiveDof<=FEM.NodeNumber);
bendingmode3=find(ActiveDof>2*FEM.NodeNumber & 3*ActiveDof<=2*FEM.NodeNumber );

ActiveBendDOF=ActiveDof(bendingmode3)-2*FEM.NodeNumber;

bendModeNo=length(bendingmode3);

%---------------Mode Shape Plot-------------------

if bendModeNo<10
    mode=bendModeNo;
else
    mode=plotmodenNo;
end

% for modeNo=  mode
    
    deformBy(ActiveBendDOF)=static_displacement(1:bendModeNo,:);
    dx=min(Xcoord):(max(Xcoord)-min(Xcoord))/100:max(Xcoord);
    dy=min(Ycoord):(max(Ycoord)-min(Ycoord))/100:max(Ycoord);
    [x3,y3]=meshgrid(dx,dy,0);
    z3=griddata(Xcoord,Ycoord,deformBy,x3,y3,'v4');
    % ZQ=griddata(Xcoord,Ycoord,Zcoord,deformUZ,x3,y3,z3,'linear');
    hf=figure;
    axes1 = axes('Parent',hf,'YTick',[0 0.8 0.801],'XTick',[0 0.6 0.601],...
        'DataAspectRatio',[1 1 1]);
    hold(axes1,'all');
    
    [max_value,max_id]=max(abs(z3(:)));
    set(gcf,'color','w')
    
    % scalefactor = 1/(max(abs(z3(:)))/max([max(abs(y3(:))) max(abs(x3(:)))]))/Stru.width;
    % 
    % if max(abs(z3(:))) == max((z3(:)))
    %     sign_max = 1;
    % elseif max(abs(z3(:))) == -min((z3(:)))
    %     sign_max = -1;
    % end
    
    fig_mode = surf(x3,y3,z3,'FaceColor','interp',...
        'LineStyle','none',...
        'EdgeColor','none',...
        'FaceLighting','phong');hold on;
    % clim([-7.145577607161163e-12 5.803395986890626e-16])
    % clim([-6.426388468125317e-11 1.869523785304594e-10])
    % clim([-4.656606151909073e-11 1.651731625297348e-12])
    % clim( [-7.178057719384481e-12 3.405127415371306e-13])
    hh=colorbar('FontSize',16)
    % colorbar('FontSize',16,'Limits',[-3e-11,0.25e-11])
    colormap(jet(20));

    if save_fig == 1
    saveas(fig_mode,[plot_folder filename_mode '_Ry'],file_typ)
    end

%%
% bendingmode=find(ActiveDof<=FEM.NodeNumber);
bendingmode4=find(ActiveDof>3*FEM.NodeNumber & ActiveDof<=4*FEM.NodeNumber );

ActiveBendDOF=ActiveDof(bendingmode4)-3*FEM.NodeNumber;

bendModeNo=length(bendingmode4);

%---------------Mode Shape Plot-------------------

if bendModeNo<10
    mode=bendModeNo;
else
    mode=plotmodenNo;
end

% for modeNo=  mode
    
    deformUX(ActiveBendDOF)=static_displacement(1:bendModeNo,:);
    dx=min(Xcoord):(max(Xcoord)-min(Xcoord))/100:max(Xcoord);
    dy=min(Ycoord):(max(Ycoord)-min(Ycoord))/100:max(Ycoord);
    [x3,y3]=meshgrid(dx,dy,0);
    z3=griddata(Xcoord,Ycoord,deformUX,x3,y3,'v4');
    % ZQ=griddata(Xcoord,Ycoord,Zcoord,deformUZ,x3,y3,z3,'linear');
    hf=figure;
    axes1 = axes('Parent',hf,'YTick',[0 0.8 0.801],'XTick',[0 0.6 0.601],...
        'DataAspectRatio',[1 1 1]);
    hold(axes1,'all');
    
    [max_value,max_id]=max(abs(z3(:)));
    set(gcf,'color','w')
    
    % scalefactor = 1/(max(abs(z3(:)))/max([max(abs(y3(:))) max(abs(x3(:)))]))/Stru.width;
    % 
    % if max(abs(z3(:))) == max((z3(:)))
    %     sign_max = 1;
    % elseif max(abs(z3(:))) == -min((z3(:)))
    %     sign_max = -1;
    % end
    
    fig_mode = surf(x3,y3,z3,'FaceColor','interp',...
        'LineStyle','none',...
        'EdgeColor','none',...
        'FaceLighting','phong');hold on;
    % clim([-7.145577607161163e-12 5.803395986890626e-16])
    % clim([-6.426388468125317e-11 1.869523785304594e-10])
    % clim([-4.656606151909073e-11 1.651731625297348e-12])
    % clim( [-7.178057719384481e-12 3.405127415371306e-13])
    hh=colorbar('FontSize',16)
    % colorbar('FontSize',16,'Limits',[-3e-11,0.25e-11])
    colormap(jet(20));

    if save_fig == 1
    saveas(fig_mode,[plot_folder filename_mode '_x'],file_typ)
    end

%%
% bendingmode=find(ActiveDof<=FEM.NodeNumber);
bendingmode5=find(ActiveDof>4*FEM.NodeNumber & ActiveDof<=5*FEM.NodeNumber );

ActiveBendDOF=ActiveDof(bendingmode5)-4*FEM.NodeNumber;

bendModeNo=length(bendingmode5);

%---------------Mode Shape Plot-------------------

if bendModeNo<10
    mode=bendModeNo;
else
    mode=plotmodenNo;
end

% for modeNo=  mode
    
    deformUY(ActiveBendDOF)=static_displacement(1:bendModeNo,:);
    dx=min(Xcoord):(max(Xcoord)-min(Xcoord))/100:max(Xcoord);
    dy=min(Ycoord):(max(Ycoord)-min(Ycoord))/100:max(Ycoord);
    [x3,y3]=meshgrid(dx,dy,0);
    z3=griddata(Xcoord,Ycoord,deformUY,x3,y3,'v4');
    % ZQ=griddata(Xcoord,Ycoord,Zcoord,deformUZ,x3,y3,z3,'linear');
    hf=figure;
    axes1 = axes('Parent',hf,'YTick',[0 0.8 0.801],'XTick',[0 0.6 0.601],...
        'DataAspectRatio',[1 1 1]);
    hold(axes1,'all');
    
    [max_value,max_id]=max(abs(z3(:)));
    set(gcf,'color','w')
    
    % scalefactor = 1/(max(abs(z3(:)))/max([max(abs(y3(:))) max(abs(x3(:)))]))/Stru.width;
    % 
    % if max(abs(z3(:))) == max((z3(:)))
    %     sign_max = 1;
    % elseif max(abs(z3(:))) == -min((z3(:)))
    %     sign_max = -1;
    % end
    
    fig_mode = surf(x3,y3,z3,'FaceColor','interp',...
        'LineStyle','none',...
        'EdgeColor','none',...
        'FaceLighting','phong');hold on;
    % clim([-7.145577607161163e-12 5.803395986890626e-16])
    % clim([-6.426388468125317e-11 1.869523785304594e-10])
    % clim([-4.656606151909073e-11 1.651731625297348e-12])
    % clim( [-7.178057719384481e-12 3.405127415371306e-13])
    hh=colorbar('FontSize',16)
    % colorbar('FontSize',16,'Limits',[-3e-11,0.25e-11])
    colormap(jet(20));

    if save_fig == 1
    saveas(fig_mode,[plot_folder filename_mode '_y'],file_typ)
    end







    %         colorbar('YLim',[-max(unique(z3)),max(unique(z3))]);
    % switch Solver
    %     case 'vibration'
    % 
    %         Natural_Freq=frequency(modeNo);
    % 
    % 
    %         %             title(['Mode shape of Mode ' num2str(modeNo) ', \omega=' num2str(Natural_Freq) 'Hz'],'FontSize',12);
    %         filename=['vibrmodeshape_mode' num2str(modeNo) ];
    % 
    % 
    % 
    %     case 'prestressed_vibr'
    % 
    %         Natural_Freq=frequency(modeNo);
    % 
    %         %             title(['Prestressed mode shape of Mode ' num2str(modeNo) ', \omega=' num2str(Natural_Freq) 'Hz'],'FontSize',12);
    % 
    %     case 'buckling'
    % 
    %         %             title(['Buckling Mode ' num2str(modeNo) ', Load Factor \lambda=' num2str(frequency(modeNo))]);%,...
    %         %             sprintf('\n'),'(\delta=' num2str(delta) ', \gamma=' num2str(gamma) ',\beta=' num2str(beta) ')'],'FontSize',12);
    % 
    %         gamma=0;
    %         beta=0;
    % 
    %         %             title(['Mode ' num2str(modeNo) ', Load Factor \lambda=' num2str(frequency(modeNo)),...
    %         %                 sprintf('\n'),'(ds/hs=' num2str(depthratio(depthNO)) ', \gamma=' num2str(gamma) ',\beta=' num2str(beta) ')'],'FontSize',12);
    % end
    %
    %
    % if exist('XXstiffener','var')
    %     for kk=1:size(XXstiffener,1)
    %         ZZstiffener(kk,:)=ones(1,size(XXstiffener(kk,:),2))*max([Stiffener.height max(abs(z3(:)))*scalefactor]);
    %         plot3(XXstiffener(kk,:),YYstiffener(kk,:),ZZstiffener(kk,:),'k-','LineWidth',3);hold on;
    %     end
    % 
    %     hold off;
    %     xlabel(['Length of the plate, a'],'FontSize',12);
    %     ylabel(['Width of the plate, b'],'FontSize',12);
    %     axis image;axis off;
    % end
    % switch Solver
    % 
    %     case 'prestressed_vibr'
    % 
    %         if lambda_b_ratio>0
    % 
    %             export_fig  ss eps
    % 
    %             filename=['vibrmodeshape_Nxy_depthraio_' num2str(depthratio(depthNO)) '_tensile'];
    % 
    %             copyfile('export_fig_out.png',[filename '.png']);
    % 
    % 
    %         elseif lambda_b_ratio<0
    % 
    %             export_fig
    % 
    %             filename=['vibrmodeshape_Nxy_depthraio_' num2str(depthratio(depthNO)) '_comp'];
    % 
    %             copyfile('export_fig_out.png',[filename '.png']);
    % 
    %         end
    % 
    % 
    % 
    %     case 'buckling'
    % 
    %         export_fig
    % 
    %         if exist('depthNO','var')
    % 
    %             filename=['buckmodeshape_Nxy_depthraio_' num2str(depthratio(depthNO)) ];
    % 
    %             copyfile('export_fig_out.png',[filename '.png']);
    %         end
    % end
    
    
    export_fig
    
%     copyfile('export_fig_out.png',[filename '.png']);
%     saveas(gcf,filename,'fig');
    %     mode_figure_name=['present_stiffenerII_mode_' num2str(plotmodenNo)];
    %     saveas(gcf,[mode_figure_name '.fig'])
    
%     saveas(fig_mode,fullfile(pic_folder, file_name_mode),'png')
%     saveas(fig_mode,fullfile(pic_folder, file_name_mode),'pdf')
%     saveas(fig_mode,fullfile(pic_folder, file_name_mode),'svg')
%     
    % % % %     %% Contour
    % % % %     %-------------------------------------------------------------
    % % % %     view(2);
    % % % %
    % % % %     %%%%% plot the contour
    % % % %     hg=figure;
    % % % %     axes2 = axes('Parent',hg,'YTick',[-1 1.8 1.801],...
    % % % %         'XTick',[-1 1.6 1.601],...
    % % % %         'DataAspectRatio',[1 1 1],...
    % % % %         'PlotBoxAspectRatio',[1 1.33333333333333 10],...
    % % % %         'LineWidth',2,...
    % % % %         'FontSize',20);
    % % % %     hold(axes2,'all');
    % % % %
    % % % %     set(gcf,'color','w')
    % % % %
    % % % %     %     [max_value,max_id]=max(abs(z3(:)));
    % % % %
    % % % %
    % % % %     contour(x3,y3,z3/z3(max_id),15,'LineWidth',3);hold on;
    % % % %     axis image;colorbar;
    % % % %
    % % % %     for kk=1:size(XXstiffener,1)
    % % % %         ZZstiffener(kk,:)=10*abs(XXstiffener(kk,:));
    % % % %         plot3(XXstiffener(kk,:),YYstiffener(kk,:),ZZstiffener(kk,:),'k-','LineWidth',3);hold on;
    % % % %     end
    % % % %     axis([0 0.601 0 0.801])
    % % % %     hold off;
    % % % %     colormap(jet);axis image;colorbar('FontSize',14);axis image; box on;
    % % % %     Z=zeros(size(FEM.elementNodes,1),size(FEM.nodesCord,2));
    % % % %     set(gcf, 'PaperPosition', [0 0 4 4]);
    % % % %     set(gcf, 'PaperSize', [4 4]);
    % % % %     %     saveas(gcf,['AIAAJ_Previbr_FreeMode_Design' Stiffener.Design '_hsbs' num2str(depthratio(depthNO))],'pdf');
% end
%%

% z3_sign = sign_max*(z3)*scalefactor;
% 
% figure
% plot(x3(1,:),z3_sign(50,:))


qr_set = [size(find(pts<=FEM.NodeNumber));
            size(find(pts>FEM.NodeNumber & pts<=2*FEM.NodeNumber))
            size(find(pts>2*FEM.NodeNumber & pts<=3*FEM.NodeNumber))
            size(find(pts>3*FEM.NodeNumber & pts<=4*FEM.NodeNumber))
            size(find(pts>4*FEM.NodeNumber & pts<=5*FEM.NodeNumber))
            ];
