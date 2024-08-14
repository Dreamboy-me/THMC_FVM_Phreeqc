function [GP,TP] = getNewPropsP(GP,GC,F,fracture)
% refer: 《Numerical Simulation of Injectivity Effects of Mineral Scaling
% and Clay Swelling in a Fractured Geothermal Reservoir》  Tianfu Xu
%   k      
% ----- =（φ/φ0）^3 * ((1-φ0 )/(1-φ))^2
%  k0   
% first can't change poro,beause using poro0 for perm

GP.rock.perm(:) = (GC.rock.poro(:)./GP.rock.poro(:)).^3.*((1-GP.rock.poro(:))./(1-GC.rock.poro(:))).^2.*GP.rock.perm(:); 
GP.rock.poro(:) =  GC.rock.poro(:);

GP.Matrix.rock.perm(:) =  GP.rock.perm(1:GP.Matrix.cells.num);
GP.Matrix.rock.poro(:) =  GP.rock.poro(1:GP.Matrix.cells.num);

Frac_cellnums =  zeros(numel(fieldnames(GP.FracGrid))+1,1);
for i = 1:numel(fieldnames(GP.FracGrid))  
    Gf = GP.FracGrid.(['Frac',num2str(i)]);  
    Frac_cellnums(i+1) = Frac_cellnums(i) + Gf.cells.num;
    Frac_cellstart = GP.Matrix.cells.num + 1 + Frac_cellnums(i);
    Frac_cellend = GP.Matrix.cells.num + Frac_cellnums(i+1);
    GP.FracGrid.(['Frac',num2str(i)]).rock.perm = GP.rock.perm(Frac_cellstart:Frac_cellend);
    GP.FracGrid.(['Frac',num2str(i)]).rock.poro = GP.rock.poro(Frac_cellstart:Frac_cellend);
end

[GP,TP] = defineNNCandTrans(GP,F,fracture);
end