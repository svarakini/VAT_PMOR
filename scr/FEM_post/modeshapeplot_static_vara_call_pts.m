
y=FEM.displacement;%    uu(:,iii);
% f=y;
%coefs = Vde(pts,:)\f(pts);
% coefs = Vde\f;
% pstar = Vde*coefs;

%     y_check=y(20020);
%     pstar_check=pstar(20020);
% 
% e1=abs((y_check-pstar_check)/y_check)*100
% e2=(norm(y-pstar))/norm(f)*100
% 
% e3=max((f-pstar)./(f+(1e-6))*100)

% plot the displacement
% modeshape_vara_qr_test
modeshapeplot_static_vara_v1
hold on

%%

NodeList=reshape(FEM.elementNodes',[size(FEM.elementNodes,1)*size(FEM.elementNodes,2),1]);
numberNodes=size(FEM.nodeCoordinates,1);

eleNodeList = NodeList+numberNodes*[0:4];


% for elem=1:FEM.numberElements
% 
%     elementDof(elem,:)=[NodeIndices NodeIndices+numberNodes NodeIndices+2*numberNodes ...
%                             NodeIndices+3*numberNodes NodeIndices+4*numberNodes]; % For the 5 DoF associated with the nodes of the element
% end                        


% load([mat_files_folder filesep 'Vde_mat_653.mat'])
[ppp]=ismember(eleNodeList,pts(1:r_stat));
[pp]=find(ppp==1);
[hh,jj]=ind2sub(size(ppp),pp);

pivot_node=unique(hh);

x_qr = FEM.nodeCoordinates_label(eleNodeList(pivot_node,1),2);
y_qr = FEM.nodeCoordinates_label(eleNodeList(pivot_node,1),3);

plot(x_qr,y_qr,'x','MarkerSize',12,'LineWidth',1,'Color','r')