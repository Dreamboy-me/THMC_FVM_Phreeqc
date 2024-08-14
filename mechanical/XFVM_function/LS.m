function [ls] = LS(iel,elem_crk)
global node element
sctr=element(iel,:);
tol_ratio=1/1000000;
x0=elem_crk(iel,1);y0=elem_crk(iel,2); 
x1=elem_crk(iel,3);y1=elem_crk(iel,4); 
for i=1:size(sctr,2)
x=node(sctr(i),1);y=node(sctr(i),2); 
phi=(y0-y1)*x+(x1-x0)*y+(x0*y1-x1*y0);
ls(i,1)=phi;
end
for j=1:size(sctr,2)
if abs(ls(j,1))<(tol_ratio*max(abs(ls)))
ls(j,1)=0;
end
end

