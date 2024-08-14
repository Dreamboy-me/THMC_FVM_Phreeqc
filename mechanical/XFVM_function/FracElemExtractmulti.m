function [Fracelemnode,elem_nodex,elem_nodey]=FracElemExtractmulti(G,F)
global node element
 
elementline=reshape(element,[],1);
elem_nodex=reshape(node(elementline,1),[],4);
elem_nodey=reshape(node(elementline,2),[],4);
% Fracture coordinates
Frac_Grid=[];
Frac_str_end=[];
for li=1:numel(F)
Frac_Grid=[Frac_Grid;G.FracGrid.(['Frac',num2str(li)]).cells.centroids];
 
nn=F(li).nodes.coords;
for i=1:size(nn,1)-1
Frac_str_end=[Frac_str_end;[nn(i,1) nn(i,2) nn(i+1,1) nn(i+1,2)]];
end

end
Fracelemnode=[];
 for j=1:size(Frac_Grid,1)
    
  for i=1:size(element,1)  
      pts=[elem_nodex(i,:)' elem_nodey(i,:)'];
      k = convhull(pts,'simplify',true);
     [in,on]=inpolygon(Frac_Grid(j,1),Frac_Grid(j,2),pts(k,1),pts(k,2));
if any(in-on)
     Fracelemnode=[Fracelemnode;[Frac_Grid(j,:),i,Frac_str_end(j,:)]]; 
end
    end
 end 
 
end