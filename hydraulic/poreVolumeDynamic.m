function [pv,sv]=poreVolumeDynamic(G)
pv_fr = poreVolumefr(G.FracGrid, G.FracGrid); % Fracture pore volume
pv_r = poreVolume(G, G.rock); % Matrix pore volume

num_fr = numel(fieldnames(G.FracGrid));
Fr_vo=cell(num_fr,1);
Fracture_cells_num=0;
for i=1:num_fr
Fracture_cells_num = Fracture_cells_num+G.FracGrid.(['Frac',num2str(i)]).cells.num;
Fracture_volumes = G.FracGrid.(['Frac',num2str(i)]).cells.volumes;
Fr_vo{i}=Fracture_volumes;
end 
Fr_pv = cell2mat(Fr_vo);
rock_cells_num=G.Matrix.cells.num; 

svfr  = Fr_pv - pv_fr; % Residual rock volume in fractures
sv    = G.cells.volumes(1:rock_cells_num) - pv_r(1:rock_cells_num) ; % Residual rock volume in matrix

pv=[pv_r(1:rock_cells_num);pv_fr];
sv=[sv;svfr];
end 