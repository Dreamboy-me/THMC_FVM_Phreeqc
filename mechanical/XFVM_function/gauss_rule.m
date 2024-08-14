function [W,Q]=gauss_rule(iel,enrich_node,elem_crk,type_elem,xTip,xVertex,normal_order,tip_order,split_order,vertex_order,junction_order,xJertex,cont)
                        
                          
global node element
sctr=element(iel,:);
tip_enr=find(enrich_node(sctr,:)==1|enrich_node(sctr,:)==11|enrich_node(sctr,:)==22);
if size(tip_enr,1)>0
 
normal_order = tip_order;
split_order  = tip_order;
vertex_order = tip_order;
junction_order=tip_order;
end


type_iel=max(type_elem(iel,:));

switch num2str(type_iel)
case {'0'} 
[W,Q]=quadrature(normal_order,'GAUSS',2);%

case{'1'} 
phi=LS(iel,elem_crk); 
 
[W,Q]=disTipQ4quad(tip_order,phi,node(sctr,:),xTip(iel,:));
 
case{'2'} 
phi=LS(iel,elem_crk); 
 
[W,Q]=discontQ4quad(split_order,phi);
 
case{'3'} 
phi=LS(iel,elem_crk); 

 
[W,Q]=disTipQ4quad(vertex_order,phi,node(sctr,:),xVertex(iel,:));
case{'4'} 
 
phi=LS(iel,elem_crk); 
if size(xJertex(iel,:),2)==8
    xJ=reshape(xJertex(iel,:),2,4);
    xJ(all(xJ== 0,2),:) = [];
   
else
xJ=reshape(xJertex(iel,:),cont,4);
xJ(all(xJ== 0,2),:) = [];

end
[W,Q]=disTipQ4quad2(junction_order,phi,node(sctr,:),xJ);
end
