function pos=posi(xCr,numnode,enrich_node)
pos=zeros(numnode,size(xCr,2));
nsnode=0;ntnode=0;njnode=0;nstnode=0;ntlnode=0;
for k=1:size(xCr,2)
    for i=1:numnode
        if(enrich_node(i,k)==2)
            pos(i,k)=(numnode+nsnode+ntnode*4+njnode+nstnode*5+ntlnode*4)+1;
            nsnode=nsnode+1;
        elseif (enrich_node(i,k)==1)
            pos(i,k)=(numnode+nsnode+ntnode*4+njnode+nstnode*5+ntlnode*4)+1;
            ntnode=ntnode+1;
        elseif(enrich_node(i,k)==3)
            pos(i,k)=(numnode+nsnode+ntnode*4+njnode+nstnode*5+ntlnode*4)+1;
            njnode=njnode+1;
        elseif(enrich_node(i,k)==22)
            pos(i,k)=(numnode+nsnode+ntnode*4+njnode+nstnode*5+ntlnode*4)+1;
            nstnode=nstnode+1;
        elseif(enrich_node(i,k)==11)
            pos(i,k)=(numnode+nsnode+ntnode*4+njnode+nstnode*5+ntlnode*4)+1;
            ntlnode=ntlnode+1;
        end
    end
end