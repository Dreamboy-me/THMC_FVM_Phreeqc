function KY=Cal_KesiYita_XFVM(Frac_elem_node,elem_nodex,elem_nodey)
global dx dy
 
KY=[];
for i=1:size(Frac_elem_node,1)
    Point(1)=Frac_elem_node(i,1);
    Point(2)=Frac_elem_node(i,2);
    
    if size(Frac_elem_node,2)==7  
        elem_num=Frac_elem_node(i,3);
        x1=elem_nodex(elem_num,1);
        x2=elem_nodex(elem_num,2);
        x3=elem_nodex(elem_num,3);
        x4=elem_nodex(elem_num,4);
        
        y1=elem_nodey(elem_num,1);
        y2=elem_nodey(elem_num,2);
        y3=elem_nodey(elem_num,3);
        y4=elem_nodey(elem_num,4);
        
    else
        x1=elem_nodex(1);
        x2=elem_nodex(2);
        x3=elem_nodex(3);
        x4=elem_nodex(4);
        
        y1=elem_nodey(1);
        y2=elem_nodey(2);
        y3=elem_nodey(3);
        y4=elem_nodey(4);
        
    end
    %µ¥Ôª×ø±ê
    
    a1   = (-x1+x2)/dx;
    a2   = (x1-x2+x3-x4)/dx/dy;
    a3   = (-x1+x4)/dy;
    a4   = x1;
    
    b1   = (-y1+y2)/dx;
    b2   = (y1-y2+y3-y4)/dx/dy;
    b3   = (-y1+y4)/dy;
    b4   = y1;
    
    %% the least squares method
    X=[0,0];
    Eps=1.0e-9;
    
    
    F =1; 
    
    %L=500;
    while(F>Eps)
        % L=L-1;
        D=0;
        F1= a4+a3*X(2)+a1*X(1)+a2*X(1)*X(2)-Point(1);
        F2= b4+b3*X(2)+b1*X(1)+b2*X(1)*X(2)-Point(2);
        F = F1*F1+F2*F2; 
        
        DF1 = a1+a2*X(2);
        DF2 = b1+b2*X(2);
        Y(1)= 2.0*(F1*DF1+F2*DF2);
        
        DF1 = a3+a2*X(1);
        DF2 = b3+b2*X(1);
        Y(2)= 2.0*(F1*DF1+F2*DF2);
        
        for j=1:2
            D=D+Y(j)*Y(j);
        end
        S=F/D;
        for k=1:2
            X(k)=X(k)-S*Y(k) ;
        end
        
    end
    Kesi(i)=X(1);
    Yita(i)=X(2);
    
    KY=[KY;[Kesi(i),Yita(i)]];
end

end 