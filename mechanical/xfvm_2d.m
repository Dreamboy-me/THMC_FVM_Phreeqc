function [u_x,u_y,nnx,nny,node] =xfvm_2d(Gmfr,L,D,nx,ny,prNew,pfNew,TrNew,TfNew,Tri,Tfi,betaT,fl,F)
%% mesh
nnx = nx+1;
nny = ny+1;
global node element dx dy
elemType = 'Q4';
pt1=[0, 0]; pt2=[L,0]; pt3=[L,D];pt4=[0,D];
dx = L/nx; dy = D/ny;
[node,element] = meshRectangularRegion(pt1, pt2, pt3, pt4, nnx,nny,elemType);
[Frac_elem_node,elem_nodex,elem_nodey]=FracElemExtractmulti(Gmfr,F);

%% define essential boundaries
uln = nnx*(nny-1)+1;       % upper left node number
urn = nnx*nny;             % upper right node number
lrn = nnx;                 % lower right node number
lln = 1;                   % lower left node number
cln = nnx*(nny-1)/2+1;     % node number at (0,0)

rightEdge= [ lrn:nnx:(uln-1); (lrn+nnx):nnx:urn ]';
leftEdge = [ uln:-nnx:(lrn+1); (uln-nnx):-nnx:1 ]';
topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';
% GET NODES ON DIRICHLET BOUNDARY AND ESSENTIAL BOUNDARY
botNodes   = unique(botEdge);     
rightNodes = unique(rightEdge);    
topNodes   = unique(topEdge);     
leftNodes  = unique(leftEdge); 
dispNodes = [leftNodes,rightNodes,botNodes,topNodes];  
dispNodes = unique(dispNodes);
ndisn = length(dispNodes);

%% PARAMETERS
nu=0.3;
E=5.0e10; 
G=E/(2*(1+nu));
lamda=nu*E/((1+nu)*(1-2*nu));
betT=betaT*E/(1-2*nu);
% Compliance matrix C
stressState='PLANE_STRAIN';
if ( strcmp(stressState,'PLANE_STRESS') )
    C = E/(1-nu^2)*[ 1   nu       0 ;
                    nu   1        0 ;
                     0   0  0.5*(1-nu) ];
else
    C = E/(1+nu)/(1-2*nu)*[ 1-nu  nu     0;
                             nu  1-nu    0;
                              0    0  0.5-nu ];
end

Cm=[1/4  1/4  0 0;
    1/4  1/4  0 0;
    1/8  3/8  0 0;
    3/8  1/8  0 0;
    0 0  1/8  3/8;
    0 0  3/8  1/8; 
    0 0  1/4  1/4;
    0 0  1/4  1/4];

Trmatrix=[-1 1 0 0;
          0 0 1 -1;
          -1 0 0 1;
          0 -1 1 0 ];
      
D=Cm*Trmatrix;
Dx_xy=D(1:4,:);%e\w\n\s
Dy_xy=D(5:8,:);%e\w\n\s

Kuu=sparse(nnx^2,nnx^2);
Kuv=sparse(nnx^2,nny^2);
Kvv=sparse(nny^2,nny^2);
Kvu=sparse(nny^2,nnx^2);

%% MATRIX EQUATIONS
 
for e=1:size(element,1)
%U方程
    e_n=element(e,:);
    uAe=[(lamda+2*G)*Dx_xy(4,:)+G*Dy_xy(2,:);-(lamda+2*G)*Dx_xy(4,:)+G*Dy_xy(1,:);...
         -(lamda+2*G)*Dx_xy(3,:)-G*Dy_xy(1,:);(lamda+2*G)*Dx_xy(3,:)-G*Dy_xy(2,:);];
    Kuu(e_n,e_n)=Kuu(e_n,e_n)+uAe;
 %UV   
    uAe=[lamda*Dy_xy(4,:)+G*Dx_xy(2,:);-lamda*Dy_xy(4,:)+G*Dx_xy(1,:);
         -lamda*Dy_xy(3,:)-G*Dx_xy(1,:);lamda*Dy_xy(3,:)-G*Dx_xy(2,:);];
    Kuv(e_n,e_n)=Kuv(e_n,e_n)+uAe; 
    
%V方程    
    uAe=[(lamda+2*G)*Dy_xy(2,:)+G*Dx_xy(4,:);(lamda+2*G)*Dy_xy(1,:)-G*Dx_xy(4,:);...
         -(lamda+2*G)*Dy_xy(1,:)-G*Dx_xy(3,:);-(lamda+2*G)*Dy_xy(2,:)+G*Dx_xy(3,:);];
    Kvv(e_n,e_n)=Kvv(e_n,e_n)+uAe; 
%VU   
    uAe=[lamda*Dx_xy(2,:)+G*Dy_xy(4,:);lamda*Dx_xy(1,:)-G*Dy_xy(4,:);...
         -lamda*Dx_xy(1,:)-G*Dy_xy(3,:);-lamda*Dx_xy(2,:)+G*Dy_xy(3,:);];
    Kvu(e_n,e_n)=Kvu(e_n,e_n)+uAe;     
 end 
 
%% fracture contribution in matrix
nf=size(Frac_elem_node,1);
Kuft=sparse(nnx^2,nf);
Kvft=sparse(nny^2,nf);

Kufn=sparse(nnx^2,nf);
Kvfn=sparse(nny^2,nf);

H=zeros(1,4);
Hp=zeros(1,4);
Hi=zeros(1,4);
Hp1=zeros(1,4);
Hi1=zeros(1,4);

KY=Cal_KesiYita_XFVM(Frac_elem_node,elem_nodex,elem_nodey);

for fr=1:nf
 
    Frac_s_e=reshape(Frac_elem_node(fr,4:end),2,2)';
    e=Frac_elem_node(fr,3); 
    e_n=element(e,:); 
    e_n_c=node(element(e,:),:); 
    seg2   = Frac_s_e(2,:)-Frac_s_e(1,:);   % tip segment
    alpha2 = atan2(seg2(2),seg2(1));  % inclination angle
   
    n=[-sin(alpha2) cos(alpha2)];%normal vector
    t=[cos(alpha2) sin(alpha2)];%shear vector 
  
     for nn=1:4
     dist=signed_distance(Frac_s_e,Gmfr.cells.centroids(e,:));
     Hp(nn)=signN(dist);
     dist=signed_distance(Frac_s_e,e_n_c(nn,:));
     Hi(nn)=signN(dist);
     H(nn)=Hp(nn)-Hi(nn);
     end 
 
      Dfx_xy=Dx_xy*H';   
      Dfy_xy=Dy_xy*H';  
% 
   
%   UUS     
      uAe1=[(lamda+2*G)*Dfx_xy(4)+G*Dfy_xy(2);-(lamda+2*G)*Dfx_xy(4)+G*Dfy_xy(1);...
          -(lamda+2*G)*Dfx_xy(3)-G*Dfy_xy(1);(lamda+2*G)*Dfx_xy(3)-G*Dfy_xy(2);]*n(1);
%   UVS
      uAe2=[lamda*Dfy_xy(4)+G*Dfx_xy(2);-lamda*Dfy_xy(4)+G*Dfx_xy(1);
          -lamda*Dfy_xy(3)-G*Dfx_xy(1);lamda*Dfy_xy(3)-G*Dfx_xy(2);]*n(2);
     Kufn(e_n,fr)=Kufn(e_n,fr)+uAe1+uAe2; 

%   VVS
      uAe1=[(lamda+2*G)*Dfy_xy(2)+G*Dfx_xy(4);(lamda+2*G)*Dfy_xy(1)-G*Dfx_xy(4);...
         -(lamda+2*G)*Dfy_xy(1)-G*Dfx_xy(3);-(lamda+2*G)*Dfy_xy(2)+G*Dfx_xy(3);]*n(2);     
%   VUS
      uAe2=[lamda*Dfx_xy(2)+G*Dfy_xy(4);lamda*Dfx_xy(1)-G*Dfy_xy(4);...
         -lamda*Dfx_xy(1)-G*Dfy_xy(3);-lamda*Dfx_xy(2)+G*Dfy_xy(3);]*n(1);
     Kvfn(e_n,fr)=Kvfn(e_n,fr)+uAe1+uAe2; 
    
  
%   UUS     
      uAe1=[(lamda+2*G)*Dfx_xy(4)+G*Dfy_xy(2);-(lamda+2*G)*Dfx_xy(4)+G*Dfy_xy(1);...
          -(lamda+2*G)*Dfx_xy(3)-G*Dfy_xy(1);(lamda+2*G)*Dfx_xy(3)-G*Dfy_xy(2);]*t(1);
%   UVS
      uAe2=[lamda*Dfy_xy(4)+G*Dfx_xy(2);-lamda*Dfy_xy(4)+G*Dfx_xy(1);
          -lamda*Dfy_xy(3)-G*Dfx_xy(1);lamda*Dfy_xy(3)-G*Dfx_xy(2);]*t(2);
      Kuft(e_n,fr)=Kuft(e_n,fr)+uAe1+uAe2; 

%   VVS
      uAe1=[(lamda+2*G)*Dfy_xy(2)+G*Dfx_xy(4);(lamda+2*G)*Dfy_xy(1)-G*Dfx_xy(4);...
         -(lamda+2*G)*Dfy_xy(1)-G*Dfx_xy(3);-(lamda+2*G)*Dfy_xy(2)+G*Dfx_xy(3);]*t(2);     
%   VUS
      uAe2=[lamda*Dfx_xy(2)+G*Dfy_xy(4);lamda*Dfx_xy(1)-G*Dfy_xy(4);...
         -lamda*Dfx_xy(1)-G*Dfy_xy(3);-lamda*Dfx_xy(2)+G*Dfy_xy(3);]*t(1);
     
     Kvft(e_n,fr)=Kvft(e_n,fr)+uAe1+uAe2; 
%}
end 

%% FRACTURE EQUATIONS
Kfft=sparse(nf,nf);
Kftn=sparse(nf,nf);
Kfnt=sparse(nf,nf);
Kffn=sparse(nf,nf);

Kfut=sparse(nf,nnx^2);
Kfvt=sparse(nf,nny^2);
Kfun=sparse(nf,nnx^2);
Kfvn=sparse(nf,nny^2);

Trfmatrix=[-1 1 0  0  0 0 0  0;
           0  0 1 -1  0 0 0  0;
           -1 0 0  1  0 0 0  0;
           0 -1 1  0  0 0 0  0;
           0  0 0  0 -1 1 0  0;
           0  0 0  0  0 0 1 -1;
           0  0 0  0 -1 0 0  1;
           0  0 0  0  0 -1 1 0;];

Lff=zeros(nf,1); 
for fr=1:nf
    fr_nodes=KY(fr,:); % Fracture coordinates
    Frac_s_e=reshape(Frac_elem_node(fr,4:end),2,2)';
    Lf=Frac_s_e(2,:)-Frac_s_e(1,:);
    Lnf=sqrt(Lf(1)^2+Lf(2)^2);%Fracture length
    Lff(fr)=Lnf;
    e=Frac_elem_node(fr,3); 
    e_n=element(e,:); 
    e_n_c=node(element(e,:),:); 
    seg2   = Frac_s_e(2,:)-Frac_s_e(1,:);   % tip segment
    alpha2 = atan2(seg2(2),seg2(1));  % inclination angle

    n=[-sin(alpha2) cos(alpha2)];%normal vector
    t=[cos(alpha2) sin(alpha2)];%shear vector
   
      DLxy_xy=[Lnf/dx*(1-fr_nodes(2)/dy) 0 0 0;
               Lnf*fr_nodes(2)/(dx*dy)   0 0 0;
               0 Lnf/dy*(1-fr_nodes(1)/dx) 0 0;
               0 Lnf*fr_nodes(1)/(dx*dy)   0 0;
               0 0 Lnf/dx*(1-fr_nodes(2)/dy) 0;
               0 0  Lnf*fr_nodes(2)/(dx*dy)  0;
               0 0 0 Lnf/dy*(1-fr_nodes(1)/dx);
               0 0 0   Lnf*fr_nodes(1)/(dx*dy);]';

      DLfxy_xy=DLxy_xy*Trfmatrix;
      DLfxy_x=DLfxy_xy(1:2,1:4);
      DLfxy_y=DLfxy_xy(3:4,5:8);
         
 
      Kfun(fr,e_n)=Kfun(fr,e_n)-lamda*DLfxy_x(1,:)-2*G*n(1)^2*DLfxy_x(1,:)-2*G*n(1)*n(2)*DLfxy_x(2,:);
      Kfvn(fr,e_n)=Kfvn(fr,e_n)-lamda*DLfxy_y(2,:)-2*G*n(2)^2*DLfxy_y(2,:)-2*G*n(1)*n(2)*DLfxy_y(1,:);
      
 
      Kfut(fr,e_n)=Kfut(fr,e_n)+G*(n(2)^2-n(1)^2)*DLfxy_x(2,:)+2*G*n(1)*n(2)*DLfxy_x(1,:);
      Kfvt(fr,e_n)=Kfvt(fr,e_n)+G*(n(2)^2-n(1)^2)*DLfxy_y(1,:)-2*G*n(1)*n(2)*DLfxy_y(2,:);
    
      for nn=1:4
      dist=signed_distance(Frac_s_e,Gmfr.cells.centroids(e,:));
      Hp(nn)=signN(dist);
      dist=signed_distance(Frac_s_e,e_n_c(nn,:));
      Hi(nn)=signN(dist);
      H(nn)=Hp(nn)-Hi(nn);
      end 
 
       Int=find(Frac_elem_node(:,3)==e);      
       if numel(Gmfr.cells.fracture.line_num{e,1})>=2&&(numel(Int)==2||numel(Int)==3||numel(Int)==4) 
       II=unique(Int);
       ff=find(Int==fr);
       Int(ff)=[];

       for i=[Int]'
       Frac_s_e_m=reshape(Frac_elem_node(i,4:end),2,2)';
       for nn=1:4        
        dist1=signed_distance(Frac_s_e_m,Gmfr.cells.centroids(e,:));
        distV1=signed_distance(Frac_s_e_m,e_n_c(nn,:));  
        Hp1(nn)=signN(dist1);
        Hi1(nn)=signN(distV1);
        
        H(nn)= H(nn)+(Hp1(nn)-Hi1(nn));  
       end      
       end    
       end   
     
      DLffxy_x=DLfxy_x*H'; 
      DLffxy_y=DLfxy_y*H';     
        
 
      Kffn(fr,fr)=Kffn(fr,fr)-(lamda*DLffxy_x(1)+2*G*n(1)^2*DLffxy_x(1)+2*G*n(1)*n(2)*DLffxy_x(2))*n(1);
      Kffn(fr,fr)=Kffn(fr,fr)-(lamda*DLffxy_y(2)+2*G*n(2)^2*DLffxy_y(2)+2*G*n(1)*n(2)*DLffxy_y(1))*n(2);      
      
      Kfnt(fr,fr)=Kfnt(fr,fr)-(lamda*DLffxy_x(1)+2*G*n(1)^2*DLffxy_x(1)+2*G*n(1)*n(2)*DLffxy_x(2))*t(1);
      Kfnt(fr,fr)=Kfnt(fr,fr)-(lamda*DLffxy_y(2)+2*G*n(2)^2*DLffxy_y(2)+2*G*n(1)*n(2)*DLffxy_y(1))*t(2);
     
      Kftn(fr,fr)=Kftn(fr,fr)+(G*(n(2)^2-n(1)^2)*DLffxy_x(2)+2*G*n(1)*n(2)*DLffxy_x(1))*n(1);
      Kftn(fr,fr)=Kftn(fr,fr)+(G*(n(2)^2-n(1)^2)*DLffxy_y(1)-2*G*n(1)*n(2)*DLffxy_y(2))*n(2);      
     
      Kfft(fr,fr)=Kfft(fr,fr)+(G*(n(2)^2-n(1)^2)*DLffxy_x(2)+2*G*n(1)*n(2)*DLffxy_x(1))*t(1);
      Kfft(fr,fr)=Kfft(fr,fr)+(G*(n(2)^2-n(1)^2)*DLffxy_y(1)-2*G*n(1)*n(2)*DLffxy_y(2))*t(2); 
end

bcwt = mean(diag(Kuu));
Kuu(dispNodes,:)= 0;
Kuu(dispNodes,dispNodes) = bcwt*speye(ndisn); 
Kuv(dispNodes,:) = 0;

bcwt = mean(diag(Kvv));
Kvv(dispNodes,:) = 0;
Kvv(dispNodes,dispNodes) = bcwt*speye(ndisn); 
Kvu(dispNodes,:) = 0;

% Stiffness matrix
K=[Kuu  Kuv  Kufn  Kuft ;
   Kvu  Kvv  Kvfn  Kvft ;
   Kfun Kfvn Kffn  Kfnt ;
   Kfut Kfvt Kftn  Kfft ;]; 
 

%% rock nodal force
f=zeros(size(K,1),1);
d = [leftNodes,rightNodes,botNodes,topNodes];
d = unique(d);
nd=1:size(node,1);
nd(d)=[];

%Equivalent nodal force
biot=1.0;
pT=biot*(prNew-0)+betT*(TrNew-Tri);

k=1;
m=0;
 for i=[nd]   
     m=m+1;
     j=i-nnx-k;
     f(i)=0.5*dy*(pT(j+1)+pT(j+1+nx)-pT(j)-pT(j+nx)); 
     f(i+nnx*nny)=0.5*dx*(pT(j+nx)+pT(j+1+nx)-pT(j+1)-pT(j)); 
     if mod(m,nnx-2)==0
     k=k+1;
     end
 end

 
pNew=reshape((pT),nx,ny);
pNew=pNew';
 
pl=0*1.0e7*ones(ny,1);
pr=0*1.0e5*ones(ny,1);
pb=pNew(1,:);
pt=pNew(end,:);

% Left boundary point
ln=leftNodes;
ln(1)=[];
ln(end)=[];
j=1;
for i=[ln]
    f(i)=(0.5*(pNew(j,1)+pNew(j+1,1))-0.5*(pl(j)+pl(j+1)))*dy;
    f(i+nnx*nny)=(0.5*(pl(j+1)+pNew(j+1,1))-0.5*(pl(j)+pNew(j,1)))*0.5*dx;
    j=j+1;
end 

%Right boundary point
rn=rightNodes;
rn(1)=[];
rn(end)=[];
j=1;
for i=[rn]
    f(i)=(0.5*(pr(j)+pr(j+1))-0.5*(pNew(j,end)+pNew(j+1,end)))*dy;
    f(i+nnx*nny)=(0.5*(pr(j+1)+pNew(j+1,end))-0.5*(pr(j)+pNew(j,end)))*0.5*dx;
    j=j+1;
end 

%Bottom boundary point
bn=botNodes;
bn(1)=[];
bn(end)=[];
j=1;
for i=[bn]
    f(i)=(0.5*(pb(j+1)+pNew(1,j+1))-0.5*(pb(j)+pNew(1,j)))*0.5*dy;
    f(i+nnx*nny)=(0.5*(pNew(1,j)+pNew(1,j+1))-0.5*(pb(j)+pb(j+1)))*dx;
    j=j+1;
end 

%Top boundary point
tn=topNodes;
tn(1)=[];
tn(end)=[];
j = 1;
for i=[tn]
    f(i)=(0.5*(pt(j+1)+pNew(end,j+1))-0.5*(pt(j)+pNew(end,j)))*0.5*dy;
    f(i+nnx*nny)=(0.5*(pt(j)+pt(j+1))-0.5*(pNew(end,j)+pNew(end,j+1)))*dx;
    j=j+1;
end 

%Top right corner point
f(nnx*nny)=(pr(end)-pNew(end,end))*0.5*dy;
 
%Bottom right corner point
f(nnx) = (pr(1)-pNew(1,end))*0.5*dy;
 
%Bottom left corner point
f(1) = -(pl(1)-pNew(1,1))*0.5*dy;

%Top left corner point
f(nnx*(nny-1)+1) = -(pl(end)-pNew(end,1))*0.5*dy;

f(dispNodes) = 0;
f(dispNodes+nnx*nny) = 0;

%% soure term

f((2*nnx*nny+1):end-nf)=(-biot*(pfNew-0)-betT*(TfNew-Tfi)).*Lff; 
 
% solve
u=K\f;
u_x = u(1:nnx*nny,1);
u_y = u(nnx*nny+1:(nnx^2+nny^2),1);

fnt = u((nnx^2+nny^2)+1:end);
fn = fnt(1:nf);
ft = fnt((nf+1):end);  
