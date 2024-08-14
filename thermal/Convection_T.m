function [Ad,dd]=Convection_T(G, velocity, thermophysical, nx, ny, cellNo, ff, cf)   
% for many fracures
Nt = G.cells.num;
Nr = G.Matrix.cells.num;
Nf = Nt-Nr;
 
vxcolumn = velocity.m(1:((nx+1)*ny));
vycolumn = velocity.m(((nx+1)*ny)+1:end);
vx = reshape(vxcolumn,nx+1,ny);
vy = reshape(vycolumn,nx,ny+1);

ix = (1+sign(vx))/2;
iy = (1+sign(vy))/2;

cprhovx = thermophysical.cpf.*vx;
cprhovy = thermophysical.cpf.*vy;
%% upwind matrix for Matrix
upw(1:nx,:) = -ix(1:nx,:).*cprhovx(1:nx,:);
ups(:,1:ny) = -iy(:,1:ny).*cprhovy(:,1:ny);
upe(1:nx,:) =  (1-ix(2:nx+1,:)) .*cprhovx(2:nx+1,:);
upn(:,1:ny) =  (1-iy(:,2:ny+1)) .*cprhovy(:,2:ny+1);

Txeast(2:nx,:)   = upe(1:nx-1,:);    Txeast(1,:)     = 0;
Tynorth(:,2:ny)  = upn(:,1:ny-1);    Tynorth(:,1)    = 0;
Txwest(1:nx-1,:) = upw(2:nx,:);      Txwest(nx,:)    = 0;
Tysouth(:,1:ny-1)= ups(:,2:ny);      Tysouth(:,ny)   = 0;

upd(1:nx,1:ny) = ix(2:nx+1,1:ny)   .*cprhovx(2:nx+1,1:ny)...
    +iy(1:nx,2:ny+1)   .*cprhovy(1:nx,2:ny+1)...
    -(1-ix(1:nx,1:ny)) .*cprhovx(1:nx,1:ny)...
    -(1-iy(1:nx,1:ny)) .*cprhovy(1:nx,1:ny);

Ds   = [Tysouth(:) Txwest(:) upd(:) Txeast(:) Tynorth(:)];

Up   = spdiags(Ds,[-nx -1 0 1 nx],nx*ny,nx*ny);

% right hand side
vb = [velocity.m;velocity.fr];
rhsr = accumarray(cellNo, -thermophysical.cpf.*vb(cf).*(ff), [Nt, 1]);

%% upwind matrix for fractures
vf =  velocity.fr;
vfr = [];
for i = 1:length(fieldnames(G.FracGrid))
    % Number of interfaces per fracture
    numface = G.FracGrid.(['Frac',num2str(i)]).faces.num;
    % Determine whether the fracture is horizontal or vertical
    if G.FracGrid.(['Frac',num2str(i)]).cells.centroids(1,1) == G.FracGrid.(['Frac',num2str(i)]).cells.centroids(2,1) % vertical
        numfry = 2*G.FracGrid.(['Frac',num2str(i)]).cells.num;
        v = vf(numfry+1:numface);
    else
        numfry = 2*G.FracGrid.(['Frac',num2str(i)]).cells.num;
        numfrx = numface-numfry;
        v = vf(1:numfrx);
    end
    
    vf = vf(numface+1:end);
    vfr = [vfr;v];
end

ifv     = (1+sign(vfr))/2;
nfrface = Nf+length(fieldnames(G.FracGrid));
Upfeast = zeros(nfrface,1);
Upfwest = zeros(nfrface,1);
updf    = zeros(nfrface,1);

A = []; B =[];

ia = 1;
ib = G.FracGrid.Frac1.cells.num;
num_fr_line = length(fieldnames(G.FracGrid));
for i = 2:num_fr_line
    A = [A ia];
    B = [B ib];
    
    ia = ib + 2;
    ib = ia+ G.FracGrid.(['Frac',num2str(i)]).cells.num -1;
    
end
A = [A ia];
B = [B ib];

for i = 1:num_fr_line
    ia = A(i);
    ib = B(i);
    
    upwf = -ifv(ia:ib).*thermophysical.cpf.*vfr(ia:ib);
    upef =  (1-ifv(ia+1:ib+1)).*thermophysical.cpf.*vfr(ia+1:ib+1);
 
    updf(ia-i+1:ib-i+1) = ifv(ia+1:ib+1).*thermophysical.cpf.*vfr(ia+1:ib+1)...
                 -(1-ifv(ia:ib)).*thermophysical.cpf.*vfr(ia:ib); 
    
    Upfwest(ia-i+1:ib-1-i+1) = upwf(2:end);
    Upfeast(ia+1-i+1:ib-i+1) = upef(1:end-1);
    
end

Dsf  = [Upfwest updf Upfeast];
Upf  = spdiags(Dsf,[-1,0,1],Nf,Nf);

method=3;
if method==1
%%  Matrix-Fracture   
nnc = G.nnc.cells(1:size(G.nnc.CI,1),:);
Vmf = velocity.mfr';
cprhoVmf = bsxfun(@times, thermophysical.cpf, Vmf);

DsVmf=cprhoVmf;
DsVmf(DsVmf<0)=0;
Upmf=cprhoVmf;
Upmf(Upmf>0)=0;

BUpmf   = sparse(nnc(:,1),nnc(:,2)-Nr,Upmf,Nr,Nf);
Bmf     = sparse(nnc(:,1),nnc(:,1),DsVmf,Nr,Nr);

%% Fracture-Matrix   

Vfm=-Vmf;
cprhoVfm = bsxfun(@times, thermophysical.cpf, Vfm);

DsVfm=cprhoVfm;
DsVfm(DsVfm<0)=0;%Upmf
Upfm=cprhoVfm;
Upfm(Upfm>0)=0;%DsVmf

BUpfm   = sparse(nnc(:,2)-Nr,nnc(:,1),-DsVfm,Nf,Nr);
Bfm     = sparse(nnc(:,2)-Nr,nnc(:,2)-Nr,DsVfm,Nf,Nf);
%% Fracture-Fracture 

nncff = G.nnc.cells(size(G.nnc.CI,1)+1:end,:);
Vff=velocity.frfr';
fnncv=[Vff;-Vff];
cprhoVff = bsxfun(@times, thermophysical.cpf, fnncv);
fnnccl=[nncff(:,1);nncff(:,2)];
fnnccr=[nncff(:,2);nncff(:,1)];
fnnc=[fnnccl,fnnccr];

DsVff=cprhoVff;
DsVff(DsVff<0)=0;
Upff=cprhoVff;
Upff(Upff>0)=0;

BUpff1   = sparse(fnnc(:,1)-Nr,fnnc(:,2)-Nr,-DsVff,Nf,Nf);
Bff      = sparse(fnnc(:,1)-Nr,fnnc(:,1)-Nr,DsVff,Nf,Nf);

BUpff1   =BUpff1+ sparse(fnnc(:,2)-Nr,fnnc(:,1)-Nr,Upff,Nf,Nf);
Bff      =Bff+ sparse(fnnc(:,2)-Nr,fnnc(:,2)-Nr,-Upff,Nf,Nf);
elseif method==2
%% Fracture-Matrix  second method 
nnc = G.nnc.cells(1:size(G.nnc.CI,1),:);
Vfm = -velocity.mfr';
cprhoVfm = bsxfun(@times, thermophysical.cpf, Vfm);

DsVfm=cprhoVfm;
DsVfm(DsVfm<0)=0;%out
Upfm=cprhoVfm;
Upfm(Upfm>0)=0;%in

%Matrix-Fracture 
Bfm    = sparse(nnc(:,2)-Nr,nnc(:,2)-Nr,-Upfm,Nf,Nf);
BUpfm  = sparse(nnc(:,2)-Nr,nnc(:,1),Upfm,Nf,Nr);
%Matrix-Fracture 
Bmf    = sparse(nnc(:,1),nnc(:,1),-Upfm,Nr,Nr);
BUpmf  = sparse(nnc(:,1),nnc(:,2)-Nr,-DsVfm,Nr,Nf);
%% Fracture-Fracture 
nncff = G.nnc.cells(size(G.nnc.CI,1)+1:end,:);
Vff=velocity.frfr';
fnncv=[Vff;-Vff];
cprhoVff = bsxfun(@times, thermophysical.cpf, fnncv);
fnnccl=[nncff(:,1);nncff(:,2)];
fnnccr=[nncff(:,2);nncff(:,1)];
fnnc=[fnnccl,fnnccr];

DsVff=cprhoVff;
DsVff(DsVff<0)=0;
Upff=cprhoVff;
Upff(Upff>0)=0;

BUpff1   = sparse(fnnc(:,1)-Nr,fnnc(:,2)-Nr,-DsVff,Nf,Nf);
Bff      = sparse(fnnc(:,1)-Nr,fnnc(:,1)-Nr,DsVff,Nf,Nf);

BUpff1   =BUpff1+ sparse(fnnc(:,2)-Nr,fnnc(:,1)-Nr,Upff,Nf,Nf);
Bff      =Bff+ sparse(fnnc(:,2)-Nr,fnnc(:,2)-Nr,-Upff,Nf,Nf);
else 
%% third method
% nnc = G.nnc.cells(1:size(G.nnc.CI,1),:);
Vmf = -velocity.fmr'; % V(matrix-->fracture) <==> -V(fracture-->matrix)
cprhoVmf = bsxfun(@times, thermophysical.cpf, Vmf); %  thermophysical.cpf-->Cpf
BUpmf    = -min(cprhoVmf,0); 
Ds       =  sum(max(cprhoVmf,0),2);
[m,n]    = size(Up);
Bmf      = spdiags(Ds,0,m,n);

Vfm      = velocity.fmr;
cprhoVfm = bsxfun(@times, thermophysical.cpf, Vfm); %  thermophysical.cpf-->Cpf
BUpfm    = -min(cprhoVfm,0);
DsT      =  sum(max(cprhoVfm,0),2); 

Vff      = velocity.vff;
cprhoVff = bsxfun(@times, thermophysical.cpf, Vff); %  thermophysical.cpf-->Cpf
BUpff1   = min(cprhoVff,0);
DsTT     = sum(max(cprhoVff,0),2);
DsT      = DsT + DsTT;
[m,n]    = size(Upf);
Bff      = spdiags(DsT,0,m,n);  % Diagonal contribution of Afm and Aff to the main diagonal of A
 
end 
%% merging matrix and fracture matrix and rhs
% Fracture-matrix solution
Upb = Up + Bmf;
Uff = Upf + Bff+BUpff1; 
Ad  = [Upb,-BUpmf;-BUpfm,Uff];
dd  = rhsr;% Concenate the RHS vectors of matrix and fracture

end 