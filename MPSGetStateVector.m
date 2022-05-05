function [stateVector] = MPSGetStateVector(matrixArray,dVector)

% Philip Daniel Blocher
% Quantum Optics Group
% Department of Physics and Astronomy, Aarhus University
% pblocher@phys.au.dk
% Algorithm for reconstructing state vector from MPS decomposition.
% Last updated: November 2017.

%% Description of input and output:
% Input:
% matrixArray returns all matrix elements making up the MPS.
% dVector is a (L times 1)-vector with the number of degrees of freedom 
% for each site.
% L is the number of sites in the 1D chain.

% Output:
% StateVector is a (n times 1)-vector with state coefficients, where n is
% the total number of states available.


%% Initialize working environment.

index = ones(size(dVector));

L = size(dVector,1);

dtot = 1;
for i = 1:1:L
    dtot = dtot * dVector(i);
end
    
stateVector = zeros(dtot,1);


%% Reconstruct state vector from MPS decomposition.

for i = 1:1:dtot
    
    % coefficient will ultimately be the state coefficient. We multiply
    % matrices according to index.    
    coefficient = matrixArray{index(1),1};
    for j = 2:1:L
        coefficient = coefficient * matrixArray{index(j),j};
    end  
   
    % Now take the trace of coefficient (which is only relevant when we
    % have periodic boundary conditions) and store the value.
    stateVector(i) = trace(coefficient);
    
    % Now update the index vector such that we calculate the next
    % coefficient.
    index(1) = index(1)+1;
    for k = 1:1:L
       if index(k) > dVector(k) && i ~= dtot
          index(k)=1;
          index(k+1)=index(k+1)+1;
       end
    end
   
end

end