function [outputMPS, stateNorm] = MPSCompressionLeft(inputMPS, dVector, dmax, eps)

% Philip Daniel Blocher
% Quantum Optics Group
% Department of Physics and Astronomy, Aarhus University
% pblocher@phys.au.dk
% Matrix Product State (MPS) compression algorithm for open boundary
% conditions (OBC).
% Last updated: November 2017.


%% Description of input and output:
% Input:
% inputMPS is a MPS, corresponding to the structure described by dVector.
% dVector is a (L times 1) vector, where each entry details the degrees of
% freedom for each site in the L-site long chain.
% dmax and eps are parameters for the compression algorithm, with dmax
% being the largest allowed matrix dimension, and eps being the cutoff
% point for weighs in the Schmidt decomposition.

% Output:
% outputMPS is the resulting MPS from the compression algorithm.
% stateNorm is the norm of the input MPS (in case this deviates from 1).

%% Initialization of environment

% A right-renormalization is made before the truncation in order to ensure 
% a normalized state.
[inputMPS,~] = MPSRenormalization(inputMPS,dVector,'RCN',[]);

% The truncation may then be performed from the left, leaving us a
% left-canonical form MPS as the result.

% The following two parameters are relevant for the compression algorithm:
L = size(dVector,1);
stateNorm = 1;

% Finally, we initialize the output MPS to be the inputMPS (such that we
% can manipulate it).
outputMPS = inputMPS; 


%% Compression algorithm
% We conduct a full-range compression sweep from left to right.
for cIndex = 1:1:L
    
    % Store all matrices at column cIndex:
    B = [];
    for rIndex = 1:1:dVector(cIndex)
        B = [B; cell2mat(outputMPS(rIndex,cIndex))];
    end
    
    % Now do SVD on B to get U0, S0, and V0Dagger.
    [U0,S0,V0] = svd(B,'econ');
    V0Dagger = V0';
    
    % In order to do the compression, we need to calculate the Schmidt
    % weights:
    S2 = S0 * S0;
    if diag(S2) == 0 % Fix for cases where S2 is the zero matrix...
        Weights = ones(size(S2));
    else
        Weights = diag(S2)./S2(1,1);
    end
    
    % We now need to select the number of dimensions to keep; this is done
    % by considering dmax and eps:
    NewDimension = sum(Weights >= eps); % might want to do a cumulative sum and make the error a sum instead.
    NewDimension = min(NewDimension,dmax);
    
    % Resize the matrices to reflect the number of dimensions to keep.
    U0 = U0(:,1:1:NewDimension);
    S0 = S0(1:1:NewDimension,1:1:NewDimension);
    V0Dagger = V0Dagger(1:1:NewDimension,:);
    
    % Calculate c0, which is the matrix to be multiplied to the right.
    c0 = S0*V0Dagger;
    
    % Now multiply c0 to the right if we are in the case cIndex ~= L or
    % store c0 as the stateNorm if we are in the case cIndex == L.
    if cIndex == L
        stateNorm = abs(c0);
        U0 = U0 * sign(c0);
    else
        for rIndex = 1:1:dVector(cIndex+1)
            outputMPS(rIndex,cIndex+1) = {c0*cell2mat(outputMPS(rIndex,cIndex+1))};
        end
    end
    
    % Finally store U0 by partitioning it into matrices with the
    % correct size. A potential problem is present here; how do we know
    % that intervalSize will be an integer? This depends on dVector(cIndex)
    % as well as sizeMeasure(2)...    
    for rIndex = 1:1:dVector(cIndex)
        sizeMeasure = size(U0);
        intervalSize = sizeMeasure(1)/dVector(cIndex);
        startIndex = (rIndex-1)*intervalSize + 1;
        endIndex = rIndex*intervalSize;
        outputMPS(rIndex,cIndex) = {U0(startIndex:1:endIndex,:)};   
    end
    
end

end

