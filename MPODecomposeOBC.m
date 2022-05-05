function [MPOArray] = MPODecomposeOBC(operatorVector,dVector)

% Philip Daniel Blocher
% Quantum Optics Group
% Department of Physics and Astronomy, Aarhus University
% pblocher@phys.au.dk
% Matrix Product Operator (MPO) decomposition algorithm for open boundary
% conditions (OBC).
% Last updated: November 2017.

%% Description of input and output:
% Input:
% operatorVector is a (n^2 times 1)-vector with state coefficients, where n is
% the total number of states available.
% dVector is a (L times 1)-vector with the number of degrees of freedom 
% for each site.
% L is the number of sites in the 1D chain.

% Output:
% MPOArray returns all matrix elements making up the MPO.

%% Initialization of working environment:

% Due to the fact that we have varying degrees of freedom on each site, we
% will create a large cell-array with dimension (L times dmax) where dmax 
% is the largest single site dimension. For computations we can then use 
% dVector to only access the appropriate number of entries.

dmax = max(dVector);
L = size(dVector,1);

% We store the matrices of the MPO for each site in the xz-plane.
MPOArray = cell(dmax, L, dmax);

% Now calculate the product of all single site dimensions; this is used
% when decomposing stateVector into matrices.
dtot = 1;
for i = 1:1:L
    di = dVector(i);
    dtot = dtot*di;
end

%% Algorithm:
c0 = operatorVector;
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
    
    xdim = rankMatrix*(di)^2;
    ydim = dtot^2;
    M0 = reshape(c0, [xdim, ydim]);
    
    % Perform SVD on the coefficient matrix M0 and calculate a new
    % coefficient vector c0:
    [U0,S0,V0] = svd(M0, 'econ');
    c0 = S0 * V0';
    
    % If we are at the end of the chain, we need to take care of the
    % operator norm:
    if i == L
        U0 = U0 .* c0;
    end
    
    % Finally, split the found U0 into matrices of relevant dimensions and
    % store them in matrixArray.    
    for j = 1:1:di
        for k = 1:1:di
            sizeMeasure = size(U0);
            intervalSizeMinor = sizeMeasure(1)/(di^2);
            intervalSizeMajor = sizeMeasure(1)/(di^1);
            startIndex = (j-1)*intervalSizeMinor + (k-1)*intervalSizeMajor + 1;
            endIndex = j*intervalSizeMinor + (k-1)*intervalSizeMajor;
            MPOArray(j,i,k) = {U0(startIndex:1:endIndex,:)};
        end
    end
   
end






end