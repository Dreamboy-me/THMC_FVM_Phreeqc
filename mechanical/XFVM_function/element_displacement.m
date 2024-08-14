function U = element_displacement(e,pos,enrich_node,u,cont)

% From the unknowns vector u, extract the parameters
% associated with the element "e"
% Then epsilon = B*U

global node element

sctr = element(e,:);
nn   = length(sctr);

% stdU contains true nodal displacement
idx = 0 ;
stdU   = zeros(2*nn,1);
for in = 1 : nn
    idx = idx + 1;
    nodeI = sctr(in) ;
    stdU(2*idx-1) = u(2*nodeI-1);
    stdU(2*idx)   = u(2*nodeI  );
end
if cont==1
    stdU=stdU;
else
    stdU=[];
end
% A contains enriched dofs

A = [];

if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    U = stdU ;
else                               % having enriched DOFs
    for in = 1 : nn
        nodeI = sctr(in) ;
        if (enrich_node(nodeI) == 2)     % H(x) enriched node
            AA = [u(2*pos(nodeI)-1);u(2*pos(nodeI))];
            A  = [A;AA];
        elseif (enrich_node(nodeI) == 1||enrich_node(nodeI)==11) % B(x) enriched node
                                          
            AA = [u(2*pos(nodeI)-1);
                u(2*pos(nodeI));
                u(2*(pos(nodeI)+1)-1);
                u(2*(pos(nodeI)+1));
                u(2*(pos(nodeI)+2)-1);
                u(2*(pos(nodeI)+2));
                u(2*(pos(nodeI)+3)-1);
                u(2*(pos(nodeI)+3));
                ];
            A  = [A;AA];
         elseif (enrich_node(nodeI) == 3)
            AA = [u(2*pos(nodeI)-1);u(2*pos(nodeI))];
            A  = [A;AA];
         elseif (enrich_node(nodeI) == 22)
              AA = [u(2*pos(nodeI)-1);
                u(2*pos(nodeI));
                u(2*(pos(nodeI)+1)-1);
                u(2*(pos(nodeI)+1));
                u(2*(pos(nodeI)+2)-1);
                u(2*(pos(nodeI)+2));
                u(2*(pos(nodeI)+3)-1);
                u(2*(pos(nodeI)+3));
                u(2*(pos(nodeI)+4)-1);
                u(2*(pos(nodeI)+4));
                ];
             A  = [A;AA]; 
        end
    end    
end

% total
U = [stdU;A];
