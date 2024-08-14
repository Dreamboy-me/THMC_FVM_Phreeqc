function [GC,TC] = getNewPropsC(GC,F,fracture,Diff)
% direct update poro and diffusion coefficient  

GC.rock.perm = GC.rock.poro * Diff;  

GC.Matrix.rock.perm(:) =  GC.rock.perm(1:GC.Matrix.cells.num);
GC.Matrix.rock.poro(:) =  GC.rock.poro(1:GC.Matrix.cells.num);

Frac_cellnums =  zeros(numel(fieldnames(GC.FracGrid))+1,1);
for i = 1:numel(fieldnames(GC.FracGrid))  
    Gf = GC.FracGrid.(['Frac',num2str(i)]);  
    Frac_cellnums(i+1) = Frac_cellnums(i) + Gf.cells.num;
    Frac_cellstart = GC.Matrix.cells.num + 1 + Frac_cellnums(i);
    Frac_cellend = GC.Matrix.cells.num + Frac_cellnums(i+1);
    GC.FracGrid.(['Frac',num2str(i)]).rock.perm = GC.rock.perm(Frac_cellstart:Frac_cellend);
    GC.FracGrid.(['Frac',num2str(i)]).rock.poro = GC.rock.poro(Frac_cellstart:Frac_cellend);
end

[GC,TC] = defineNNCandTrans(GC,F,fracture);
end