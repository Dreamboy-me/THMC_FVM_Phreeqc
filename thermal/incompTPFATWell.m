function state = incompTPFATWell(state, G, T, fluidT, fluidP, velocity, prNew, TrOld, TfOld, rhof, muf, thermophysical, Cpm, Cpf, pv, sv, nx, ny, dt, Tinj, varargin)                                            
% Solve heat transferproblem using TPFA method.
 
%   This function uses a two-point flux approximation (TPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.
 
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity,
%   then the input values 'xr' and 'xw' are returned unchanged and a
%   warning is printed in the command window. This warning is printed with
%   message ID 

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of  the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


   opt = struct('bc', [], 'src', [], 'wells', [], ...
                'Wells', [], 'bcp',[],...
                'LinSolve',     @mldivide,        ...
                'MatrixOutput', false, ...
                'Verbose',      mrstVerbose,...
                'gravity',      gravity(), ...
                'condition_number',false,...
                'pc_form','nonwetting',...
                'reduce',false,...
                'use_trans',false);

   opt = merge_options(opt, varargin{:});
 
   g_vec   = opt.gravity;
   % If gravity is overriden, we cannot say anything about the effects of
   % gravity on rhs.
   grav = (norm(g_vec(1:G.griddim)) > 0) || isfield(G, 'grav_pressure');

   if all([~opt.MatrixOutput , ...
           isempty(opt.bc)   , ...
           isempty(opt.src)  , ...
           isempty(opt.bcp)  , ...
           isempty(opt.Wells), ~grav]) % original:isempty(opt.wells)
      warning(msgid('DrivingForce:Missing'),                   ...
             ['No external driving forces present in model--', ...
              'state remains unchanged.\n']);
   end

   % ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Setting up linear system...\t\t\t');
   t0 = ticif (opt.Verbose);

   % Preliminaries
   [neighborship, n_isnnc] = getNeighbourship(G, 'Topological', true);
   [cellNo, cf, cn_isnnc] = getCellNoFaces(G);
   nif    = size(neighborship, 1);
   ncf    = size(cf, 1);
   nc     = G.cells.num;
   nw     = length(opt.wells);
   n      = nc + nw;

   [mobr, mobf, omegaf, rhor] = dynamic_quantities(state, fluidT, rhof, muf);
   totmob = sum(mobr, 2);
   
   % Compute effective (mobility-weighted) transmissibilities.
   [T, ft] = compute_trans(G, T, cellNo, cf, neighborship, totmob, opt);

   % Identify internal faces
   i  = all(neighborship ~= 0, 2);

   % Boundary conditions and source terms.
   hh = zeros(nif, 1);
   dF = false(nif, 1);
   [grav, ff] = deal(zeros(ncf, 1));

   [ff(~cn_isnnc), gg, hh(~n_isnnc), ...
    grav(~cn_isnnc), dF(~n_isnnc), dC] = ...
      computePressureRHS(G, omegaf, opt.bc, opt.src);

   % made to add capillary pressure
   pc = getIncompCapillaryPressure(state, fluidP);
   if ~isempty(pc)
      if any(abs(pc) > 0)
         cc = capPressureRHS(G, mobf, pc, opt.pc_form);
         grav = grav + cc;
      end
   end

   sgn = 2*(neighborship(cf, 1) == cellNo) - 1;
   j   = i(cf) | dF(cf);
   fg  = accumarray(cf(j), grav(j).*sgn(j), [nif, 1]);
   if ~isempty(opt.bcp)
      fg(opt.bcp.face) = fg(opt.bcp.face) + opt.bcp.value;
      warning('mrst:periodic_bc', ...
             ['Face pressures are not well defined for ', ...
              'periodic boundary faces']);

      if any(G.faces.neighbors(:,1) == G.faces.neighbors(:,2))
         error(['Periodic boundary: This code do not work ', ...
                'if a face is in and outflow']);
      end
   end
 
   rhocpeff_m = pv(1:G.Matrix.cells.num).*rhof(1:G.Matrix.cells.num).* Cpf + sv(1:G.Matrix.cells.num).*rhor.* Cpm;
   rhocpeff_fr = pv(G.Matrix.cells.num+1:end).*rhof(G.Matrix.cells.num+1:end).* Cpf + sv(G.Matrix.cells.num+1:end).*rhor.* Cpm;
 
   rhocpeff     = [rhocpeff_m; rhocpeff_fr];
   
   [Ad,dd] = Convection_T(G, velocity, thermophysical, nx, ny, cellNo(~cn_isnnc), ff(~cn_isnnc), cf(~cn_isnnc));  

   % Wells --------------------------
  
   if~isempty(opt.Wells)
       
       w =opt.Wells;
       ww=zeros(nc,1);
       is_press = strcmpi('bhp', w(1).type);
       if  is_press
           qinj=w(1).WI*rhof(w(1).cells)/muf(w(1).cells)*(prNew(w(1).cells)-w(1).val); 
           
           if qinj<0
               qinj=-qinj;
           end
           ww(w(1).cells)=thermophysical.cpf*Tinj*qinj;
       else
           ww(w(1).cells)=rhof(w(1).cells).*thermophysical.cpf*Tinj*w(1).val;
       end
       
       %production well
       wq=zeros(nc,1);
       is_press = strcmpi('rate', w(2).type);
       if is_press
           wq(w(2).cells)=rhof(w(2).cells)*thermophysical.cpf*w(2).val;
       else
           qpro=w(2).WI*rhof(w(2).cells)/muf(w(2).cells)*(prNew(w(2).cells)-w(2).val);%production flux
           
           wq(w(2).cells)=thermophysical.cpf*qpro;
       end
   else
       ww=zeros(nc,1);
       wq=zeros(nc,1);
       
   end
  
   %-----------------------------------------
   
   TOld = [TrOld;TfOld];
   rhs = dd + ...
         accumarray(cellNo, -ft(cf).*(sgn.*fg(cf)+ff), [n, 1]) + ...
         [gg; zeros(nw, 1)]                                    + ...
         accumarray(cellNo, -hh(cf), [n, 1])                   + ...
         [rhocpeff.*TOld./dt; zeros(nw, 1)] + ...
         ww; clear sgn;  
 

   % Add up internal face transmissibilities plus Dirichlet pressure
   % faces for each cell.
   d  = zeros(G.cells.num, 1);   
   d  = d + ...
       accumarray(cellNo(dF(cf)), T(dF(cf)), [nc, 1]) +...
       accumarray(reshape(neighborship(i,:), [], 1), ...
       repmat(ft(i), [2,1]),  [nc, 1])                +...
       wq + ...
       ww./Tinj + ...
       [rhocpeff./dt; zeros(nw, 1)];  
   
   % Assemble coefficient matrix for internal faces.  Boundary conditions
   % may introduce additional diagonal entries.  Also, wells introduce
   % additional equations and unknowns.
   I  = [neighborship(i,1); neighborship(i,2); (1:nc)'];
   J  = [neighborship(i,2); neighborship(i,1); (1:nc)'];
   V  = [-ft(i); -ft(i); d]; clear d;
    
   A  = sparse(double(I), double(J), V, nc, nc);  
   A = A + Ad; 
 
   if(~opt.reduce)
    clear I J V C D;
   end
   tocif(opt.Verbose, t0);
   
   % if reduce to cells


   % ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Solving linear system...\t\t\t');
   t0 = ticif (opt.Verbose);
 
   if opt.condition_number
      disp('***************************************');
      disp(['Conditon number is :   ', num2str(condest(A))]);
      disp('***************************************');
   end
   if opt.MatrixOutput
      state.A   = A;
      state.rhs = rhs;
   end
   
   if(opt.reduce)
      rhs_r = rhs(1:nc);
      A_r = A(1:nc,1:nc);
      rhs_w = rhs(nc+1:end);
      rhs_r = rhs_r - C'*(D\ rhs_w);
      A_r = A_r - C'*(D\C);
      state.A = A_r;
      state.rhs = rhs_r;
      pT_r = opt.LinSolve(A_r, rhs_r);
      pT = nan(size(rhs));
      pT(1:nc) = pT_r;
      pT(nc+1:end) = D\(rhs_w - C*pT_r);
   else 
    pT = opt.LinSolve(A, rhs);
   end 
   tocif(opt.Verbose, t0);

   clear A rhs;

   % ---------------------------------------------------------------------
   dispif(opt.Verbose, 'Computing fluxes, face pressures etc...\t\t');
   t0 = ticif (opt.Verbose);

   % Reconstruct face pressures and fluxes.
   ftemperature     =  ...
          accumarray(cf, (pT(cellNo)+grav).*T, [nif, 1])./ ...
          accumarray(cf(:,1), T, [nif,1]);


   % Neumann faces
   b         = any(G.faces.neighbors==0, 2);
   ftemperature(b) = ftemperature(b) - hh(b)./ft(b);  

   % Dirichlet faces
   ftemperature(dF) = dC;


   % Sign for boundary faces
   sgn  = 2*(G.faces.neighbors(~i,2)==0)-1;
   ni   = neighborship(i,:);
   fluxT = -accumarray(find(i),  ft(i) .*(pT(ni(:,2))-pT(ni(:,1))-fg(i)), [nif, 1]);
   c    = sum(G.faces.neighbors(~i,:),2) ;
   fg  = accumarray(cf, grav, [nif, 1]);
   fluxT(~i) = -sgn.*ft(~i).*( ftemperature(~i) - pT(c) - fg(~i) );
 
   state.temperature = pT(1 : nc); 
   state.TrNew = pT(1 : G.Matrix.cells.num);  
   state.TfNew = pT(G.Matrix.cells.num+1:nc); 
   state.fluxT         = fluxT;
   state.facetemperature     = ftemperature;   
   
   if ~isempty(opt.Wells)
       W = opt.Wells;
       nw = length(opt.Wells);
       for k = 1:nw
           state.wellSol(k) = struct('thermalpower',0,'temperature',0);
           if W(k).Type_well == "injection"
               state.wellSol(k).thermalpower = qinj * thermophysical.cpf * Tinj;
               state.wellSol(k).temperature = Tinj;
           else
               state.wellSol(k).thermalpower = qpro * thermophysical.cpf * pT(w(k).cells);
               state.wellSol(k).temperature = pT(w(k).cells);
           end
       end
   end

   tocif(opt.Verbose, t0);
end

%--------------------------------------------------------------------------

function [mobr, mobf, omegaf, rhor] = dynamic_quantities(state, fluid, rhof, muf)
   [rhor, kr, mur] = getIncompProps(state, fluid);
   mobr    = bsxfun(@rdivide, kr, mur);
   mobf    = bsxfun(@rdivide, kr, muf);
   totmobr = sum(mobr, 2);
   omegaf  = sum(bsxfun(@times, mobr, rhof), 2) ./ totmobr;
end

%--------------------------------------------------------------------------

function [T, ft] = compute_trans(G, T, cellNo, cellFaces, neighborship, totmob, opt) %#ok
    niface = size(neighborship, 1);
    if opt.use_trans
      neighborcount = sum(neighborship > 0, 2);
      assert (numel(T) == niface, ...
             ['Expected one transmissibility for each interface ', ...
              '(=%d) but got %d.'], niface, numel(T));

      fmob = accumarray(cellFaces, totmob(cellNo), ...
                        [niface, 1]);
  
      fmob = fmob ./ neighborcount;
      ft   = T .* fmob;

      % Synthetic one-sided transmissibilities.
      th = ft .* neighborcount;
      T  = th(cellFaces(:,1));

   else

      % Define face transmissibility as harmonic average of mobility
      % weighted one-sided transmissibilities.
      %
      assert (numel(T) == numel(cellNo), ...
             ['Expected one one-sided transmissibility for each ', ...
              'half face (=%d), but got %d.'], numel(cellNo), numel(T));

      T  = T .* totmob(cellNo);
      ft = 1 ./ accumarray(cellFaces, 1 ./ T, [niface, 1]);

   end
end
