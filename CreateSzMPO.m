function [Sz_MPO] = CreateSzMPO(dVector,index)
% Philip Daniel Blocher
% Quantum Optics Group
% Department of Physics and Astronomy, Aarhus University
% pblocher@phys.au.dk
% Last updated: January 2019.

%% Creation of a Base MPO with delta functions in every entry
L = size(dVector,1);
dmax = max(dVector);
BaseMPO = cell(dmax,L,dmax);

for l = 1:1:L
    dl = dVector(l);
    
    for i = 1:1:dl
        for j = 1:1:dl
            BaseMPO{i,l,j} = (i==j); % Assigns a delta function to each entry.
        end
    end
end


%% Create the Sz operator on the index site 
dindex = dVector(index);

Sz_MPO = BaseMPO;

for i = 1:1:dindex
    for j = 1:1:dindex        
        Sz_MPO{i,index,j} = 0;
    end
end

Sz_MPO{1,index,1} = 1;
Sz_MPO{2,index,2} = -1;

end