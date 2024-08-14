function [ival] = point_in_line(xTip,q1,q2)
 
AC=xTip(1,:)-q1(1,:);
BC=xTip(1,:)-q2(1,:);
AB=q2(1,:)-q1(1,:);
ab=sqrt(AB(1,1)^2+AB(1,2)^2);
ac=sqrt(AC(1,1)^2+AC(1,2)^2);
bc=sqrt(BC(1,1)^2+BC(1,2)^2);
r=AC(1,1)*AB(1,2)-AC(1,2)*AB(1,1);
if r==0&&((ac+bc)-ab<1.0e-8)
    ival=1;
else 
    ival=0;
end
    
end

