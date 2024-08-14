function [GT,TT] = getNewPropsT(GT,GC,F,fracture,cond_fluid,cond_solid)
% different from GP, first change poro

GT.rock.poro(:) =  GC.rock.poro(:);
GT.rock.perm = GT.rock.poro * cond_fluid + (1-GT.rock.poro) * cond_solid;

GT.Matrix.rock.perm(:) =  GT.rock.perm(1:GT.Matrix.cells.num);
GT.Matrix.rock.poro(:) =  GT.rock.poro(1:GT.Matrix.cells.num);

Frac_cellnums =  zeros(numel(fieldnames(GT.FracGrid))+1,1);
for i = 1:numel(fieldnames(GT.FracGrid))  
    Gf = GT.FracGrid.(['Frac',num2str(i)]);  
    Frac_cellnums(i+1) = Frac_cellnums(i) + Gf.cells.num;
    Frac_cellstart = GT.Matrix.cells.num + 1 + Frac_cellnums(i);
    Frac_cellend = GT.Matrix.cells.num + Frac_cellnums(i+1);
    GT.FracGrid.(['Frac',num2str(i)]).rock.perm = GT.rock.perm(Frac_cellstart:Frac_cellend);
    GT.FracGrid.(['Frac',num2str(i)]).rock.poro = GT.rock.poro(Frac_cellstart:Frac_cellend);
end

[GT,TT] = defineNNCandTrans(GT,F,fracture);
end