function [outputMPS, stateNorm] = MPSRenormalization(inputMPS,dVector,Mode,TargetColumn)

% Philip Daniel Blocher
% Quantum Optics Group
% Department of Physics and Astronomy, Aarhus University
% pblocher@phys.au.dk

% Algorithm for renormalization of matrix product states (MPS).
% Last updated: November 2017.


%% Description of input and output:
% Input:
% inputMPS is a matrix product state.
% dVector is a (L times 1)-vector with the number of degrees of freedom 
% for each site.
% L is the number of sites in the 1D chain.

% Mode = 'RCN', 'LCN', 'MixedL', 'MixedR'
% If 'MixedL','MixedR' chosen, TargetColumn must also be declared.
% 'MixedL' renormalizes the matrices 1 through TargetColumn
% 'MixedR' renormalizes the matrices TargetColumn through end.
% 'LCN' performs a left-canonical renormalization.
% 'RCN' performs a right-canonical renormalization.

% Output:
% outputMPS is the renormalized version of inputMPS.
% stateNorm returns the norm of the input state; outputMPS is renormalized.


%% Mode Selector
% The following logical statement allows the user to call a specific
% version of the renormalization algorithm.

% Length of chain
L = size(dVector,1);

if strcmp(Mode,'RCN')
    [outputMPS, stateNorm] = RightNormalization(inputMPS,dVector,1);
elseif strcmp(Mode, 'LCN')
    [outputMPS, stateNorm] = LeftNormalization(inputMPS,dVector,L);
elseif strcmp(Mode, 'MixedL')
    [outputMPS, stateNorm] = LeftNormalization(inputMPS,dVector,TargetColumn);
elseif strcmp(Mode, 'MixedR')
    [outputMPS, stateNorm] = RightNormalization(inputMPS,dVector,TargetColumn);
end

end


function [matrixStorage, stateNorm] = RightNormalization(inputMatrixStorage,dVector,TargetColumn)
%% Initialization of Environment
L = size(dVector,1); % Length of chain
matrixStorage = inputMatrixStorage; % Initialize output to structure identical to input
stateNorm = 0; % This will change if relevant.

%% Full-range sweep from the right
for cIndex = L:-1:TargetColumn
    B = [];
    
    for rIndex = 1:1:dVector(cIndex)
        B = [B, cell2mat(matrixStorage(rIndex,cIndex))];
    end
    
    [U0,S0,V0] = svd(B,'econ');
    V0Dagger = V0';
    c0 = U0*S0;
    
    if cIndex == 1 % This particular if-condition fixes wrongful results in case of lack of normalized state vector.
        stateNorm = abs(c0);
        V0Dagger = V0Dagger*sign(c0);
    else % Multiply U0*S0 to the left
        for rIndex = 1:1:dVector(cIndex-1)
            matrixStorage(rIndex,cIndex-1) = {cell2mat(matrixStorage(rIndex,cIndex-1))*c0};
        end
       
    end
    
    for rIndex = 1:1:dVector(cIndex)
        sizeMeasure = size(V0Dagger);
        intervalSize = sizeMeasure(2)/dVector(cIndex);
        startIndex = (rIndex-1)*intervalSize + 1;
        endIndex = rIndex*intervalSize;
        matrixStorage(rIndex,cIndex) = {V0Dagger(:,startIndex:1:endIndex)};   
    end
    
end

end

function [matrixStorage,stateNorm] = LeftNormalization(inputMatrixStorage,dVector,TargetColumn)
%% Initialization of Environment
L = size(dVector,1); % Length of chain
matrixStorage = inputMatrixStorage; % Initialize output to structure identical to input
stateNorm = 0; % This will change if relevant.

%% Normalization sweep from the left to TargetColumn
for cIndex = 1:1:TargetColumn
    A = [];
    
    for r = 1:1:dVector(cIndex)
        A = [A; cell2mat(matrixStorage(r,cIndex))];
    end
    
    [U0,S0,V0] = svd(A, 'econ');
    c0 = S0*(V0');
    
    if cIndex == L % This particular if-condition fixes wrongful results in case of lack of normalized state vector.
        stateNorm = abs(c0);
        U0 = U0 * sign(c0);
    else % Multiply S0*V0' to the right
        for r = 1:1:dVector(cIndex+1)
            matrixStorage(r,cIndex+1) = {c0*cell2mat(matrixStorage(r,cIndex+1))};
        end
       
    end
    
    for r = 1:1:dVector(cIndex)
        sizeMeasure = size(U0);
        intervalSize = sizeMeasure(1)/dVector(cIndex);
        startIndex = (r-1)*intervalSize + 1;
        endIndex = r*intervalSize;
        matrixStorage(r,cIndex) = {U0(startIndex:1:endIndex,:)};   
    end
    
end

end