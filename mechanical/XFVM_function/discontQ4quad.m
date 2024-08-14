function [W,Q]=discontQ4quad(order,phi)
% [W,Q]=discontQ4quad(order,phi,elem_crk,nodes)
 

corner = [1 2 3 4 1];
node   = [-1 -1; 1 -1; 1 1; -1 1];

% ntip=G_L_coor(elem_crk(1),elem_crk(2),nodes);
% ntip1=G_L_coor(elem_crk(3),elem_crk(4),nodes);
%  node=[node;ntip;ntip1];
%  nodes=[nodes;elem_crk(1),elem_crk(2);elem_crk(3),elem_crk(4)];
%  tri = delaunay(nodes(:,1),nodes(:,2))
%  triplot(tri,nodes(:,1),nodes(:,2));
% pause
% loop on element edges
for i = 1 : 4
    n1 = corner(i);
    n2 = corner(i+1);
    if ( phi(n1)*phi(n2) < 0 )
        r    = phi(n1)/(phi(n1)-phi(n2));
        pnt  = (1-r)*node(n1,:)+r*node(n2,:); 
        node = [node;pnt];
    end
    if(abs(phi(n1)*phi(n2) )< 1.0e-4)
        pnt=node(n1,:);
        node=[node;pnt];
    end
    
end


if size(find(phi==0),1)==2 
    node=[node;0 0];
end
% get decompused triangles
tri = delaunay(node(:,1),node(:,2)); 
tri = tricheck(node,tri);


 
% loop over subtriangles to get quadrature points and weights
pt = 1;
for e = 1:size(tri,1)
    [w,q]=quadrature(order,'TRIANGULAR',2);
    % transform quadrature points into the parent element
     
    coord = node(tri(e,:),:);
    a = det([coord,[1;1;1]])/2;
    if ( a<0 )  % need to swap connectivity
        coord = [coord(2,:);coord(1,:);coord(3,:)];
        a = det([coord,[1;1;1]])/2;
    end

   if ( a~=0 )
%   if (a>1e-8)
        for n=1:length(w)
            N=lagrange_basis('T3',q(n,:));
            Q(pt,:) = N'*coord;
            W(pt,1) = 2*w(n)*a;
            pt = pt+1;
        end
   end

end








