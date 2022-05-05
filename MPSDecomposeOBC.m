function [MPSArray, stateNorm] = MPSDecomposeOBC(stateVector,dVector)

% Philip Daniel Blocher
% Quantum Optics Group
% Department of Physics and Astronomy, Aarhus University
% pblocher@phys.au.dk
% Matrix Product State (MPS) decomposition algorithm for open boundary
% conditions (OBC).
% Last updated: November 2017.

%% Description of input and output:
% Input:
% stateVector is a (n times 1)-vector with state coefficients, where n is
% the total number of states available.
% dVector is a (L times 1)-vector with the number of degrees of freedom 
% for each site.

% Output:
% MPSArray returns all matrix elements making up the MPS.
% stateNorm returns the norm of the input state; MPS is normalized.

%% Initialization of working environment:

% Due to the fact that we have varying degrees of freedom on each site, we
% will create a large cell-array with dimension (L times dmax) where dmax 
% is the largest single site dimension. For computations we can then use 
% dVector to only access the appropriate number of entries.

dmax = max(dVector);
L = size(dVector,1);

MPSArray = cell(dmax, L, 1);

% Now calculate the product of all single site dimensions; this is used
% when decomposing stateVector into matrices.
dtot = 1;
for i = 1:1:L
    di = dVector(i);
    dtot = dtot*di;
end


%% Algorithm:
c0 = stateVector;
rankMatrix = 1;

for i = 1:1:L
    
    % Load the i'th site's degrees of freedom and update dtot to represent
    % the remaining (by division).    
    di = dVector(i);
    dtot = dtot / di;
    
    % Calculate xdim and ydim, the dimensions that our coefficient vector
    % should be reshapen into.
    if i ~= 1
        rankMatrix = size(U0,2);
    end
    
    xdim = rankMatrix*di;
    ydim = dtot;
    
    
    M0 = reshape(c0, [xdim, ydim]);
    
    % Perform SVD on the coefficient matrix M0 and calculate a new
    % coefficient vector c0:
    [U0,S0,V0] = svd(M0, 'econ');
    c0 = S0 * V0';
    
    % If we are at the end of the chain, we need to store the stateNorm and
    % return it.
    if i == L
        stateNorm = abs(c0);
        U0 = U0 .* sign(c0);
    end
    
    % Finally, split the found U0 into matrices of relevant dimensions and
    % store them in matrixArray.    
    for j = 1:1:di
        sizeMeasure = size(U0);
        intervalSize = sizeMeasure(1)/di;
        startIndex = (j-1)*intervalSize + 1;
        endIndex = j*intervalSize;
        MPSArray(j,i) = {U0(startIndex:1:endIndex,:)};     
    end
   
end

end