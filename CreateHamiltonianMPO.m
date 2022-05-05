function [MPO_Odd, MPO_Even] = CreateHamiltonianMPO(hB,hLE,hRE,tau,d,L)

% Philip Daniel Blocher #201205991
% September 2015
% Department of Physics, Aarhus University
% philipblocher@gmail.com
% Creates an MPO for a given Hamiltonian for use in odd-even time
% evolution schemes.

%% Creation of the Odd and Even bond MPOs

% Be advised that hB denotes the bulk Hamiltonian, hLE denotes the Left
% Edge Hamiltonian, and hRE denotes the Right Edge Hamiltonian. In the
% following, B will denote bulk, LE will denote Left Edge, and RE will
% denote Right Edge.

% Matrix Exponentiation
UB = expm(-1i*hB*tau);
ULE = expm(-1i*hLE*tau);
URE = expm(-1i*hRE*tau);

% The regrouping of indices of the time evolution operators will now be
% done through a number of for loops.
MB = NaN(size(UB));
MLE = MB;
MRE = MB;

for sigma1 = 1:1:d;
    for sigma1p = 1:1:d;
        for sigma2 = 1:1:d;
            for sigma2p = 1:1:d
                alphaBefore = (sigma1 - 1)*d + sigma2;
                betaBefore = (sigma1p - 1)*d + sigma2p;
                alphaAfter = (sigma1 - 1)*d + sigma1p;
                betaAfter = (sigma2 - 1)*d + sigma2p;
                
                MB(alphaAfter,betaAfter) = UB(alphaBefore,betaBefore);
                MLE(alphaAfter,betaAfter) = ULE(alphaBefore,betaBefore);
                MRE(alphaAfter,betaAfter) = URE(alphaBefore,betaBefore);
            end
        end
    end
end

clear('UB')
clear('URE')
clear('ULE')
% SVDs are now carried out on the regrouped matrices MB, MLE, and MRE.
[WB, SB, VB] = svd(MB,'econ');
UB = WB*sqrt(SB);
UBbar = sqrt(SB) * ctranspose(VB);
clear('WB')
clear('SB')
clear('VB')

[WLE, SLE, VLE] = svd(MLE,'econ');
ULE = WLE * sqrt(SLE);
ULEbar = sqrt(SLE) * ctranspose(VLE);
clear('WLE')
clear('SLE')
clear('VLE')

[WRE, SRE, VRE] = svd(MRE,'econ');
URE = WRE * sqrt(SRE);
UREbar = sqrt(SRE) * ctranspose(VRE);
clear('WRE')
clear('SRE')
clear('VRE')

% The found matrices are split into row- and column vectors
MatrixStorageB = cell(d,2,d);
MatrixStorageLE = cell(d,2,d);
MatrixStorageRE = cell(d,2,d);
for sigma = 1:1:d;
    for sigmaprime = 1:1:d;
        MatrixStorageB{sigma,1,sigmaprime} = UB((sigma-1)*d + sigmaprime,:);
        MatrixStorageB{sigma,2,sigmaprime} = UBbar(:,(sigma-1)*d + sigmaprime);
        
        MatrixStorageLE{sigma,1,sigmaprime} = ULE((sigma-1)*d + sigmaprime,:);
        MatrixStorageLE{sigma,2,sigmaprime} = ULEbar(:,(sigma-1)*d + sigmaprime);
        
        MatrixStorageRE{sigma,1,sigmaprime} = URE((sigma-1)*d + sigmaprime,:);
        MatrixStorageRE{sigma,2,sigmaprime} = UREbar(:,(sigma-1)*d + sigmaprime);
    end
end

clear('UB'); clear('UBbar'); clear('ULE'); clear('ULEbar');
clear('URE'); clear('UREbar');

% The following code builds the MPOs - discrimination is done to whether L
% is even or odd.

% Initialization of final storage
MPO_Even = cell(d,L,d);
MPO_Odd = cell(d,L,d);

site = 1;
for sigma = 1:1:d;
    for sigmaprime = 1:1:d;
        MPO_Odd{sigma,site,sigmaprime} = MatrixStorageLE{sigma,1*(mod(site,2)==1)+2*(mod(site,2)==0),sigmaprime}; % Left edge Hamiltonian
        MPO_Even{sigma,site,sigmaprime} = (sigma==sigmaprime); % Identity matrix
    end
end

if L > 2
    site = 2;
    for sigma = 1:1:d;
        for sigmaprime = 1:1:d;
            MPO_Odd{sigma,site,sigmaprime} = MatrixStorageLE{sigma,1*(mod(site,2)==1)+2*(mod(site,2)==0),sigmaprime}; % Left edge Hamiltonian
            MPO_Even{sigma,site,sigmaprime} = MatrixStorageB{sigma,2*(mod(site,2)==1)+1*(mod(site,2)==0),sigmaprime}; % Even bulk Hamiltonian
        end
    end
    
    if L > 3
        for site = 3:1:(L-2);
            for sigma = 1:1:d;
                for sigmaprime = 1:1:d;
                    MPO_Odd{sigma,site,sigmaprime} = MatrixStorageB{sigma,1*(mod(site,2)==1)+2*(mod(site,2)==0),sigmaprime}; % Odd bulk Hamiltonian
                    MPO_Even{sigma,site,sigmaprime} = MatrixStorageB{sigma,2*(mod(site,2)==1)+1*(mod(site,2)==0),sigmaprime}; % Even bulk Hamiltonian
                end
            end
        end
    end

    site = L-1;
    if mod(L,2) == 0
        for sigma = 1:1:d;
            for sigmaprime = 1:1:d;
                MPO_Odd{sigma,site,sigmaprime} = MatrixStorageRE{sigma,1*(mod(site,2)==1)+2*(mod(site,2)==0),sigmaprime}; % Right edge Hamiltonian
                MPO_Even{sigma,site,sigmaprime} = MatrixStorageB{sigma,2*(mod(site,2)==1)+1*(mod(site,2)==0),sigmaprime}; % Even bulk Hamiltonian
            end
        end
    elseif mod(L,2) == 1
         for sigma = 1:1:d;
            for sigmaprime = 1:1:d;
                MPO_Odd{sigma,site,sigmaprime} = MatrixStorageB{sigma,1*(mod(site,2)==1)+2*(mod(site,2)==0),sigmaprime}; % Odd bulk Hamiltonian
                MPO_Even{sigma,site,sigmaprime} = MatrixStorageRE{sigma,2*(mod(site,2)==1)+1*(mod(site,2)==0),sigmaprime}; % Right edge Hamiltonian
            end
        end   
    end

end

site = L;
if mod(L,2) == 0
    for sigma = 1:1:d;
        for sigmaprime = 1:1:d;
            MPO_Odd{sigma,site,sigmaprime} = MatrixStorageRE{sigma,1*(mod(site,2)==1)+2*(mod(site,2)==0),sigmaprime}; % Right edge Hamiltonian
            MPO_Even{sigma,site,sigmaprime} = (sigma==sigmaprime); % Identity matrix
        end
    end
elseif mod(L,2) == 1
     for sigma = 1:1:d;
        for sigmaprime = 1:1:d;
            MPO_Odd{sigma,site,sigmaprime} = (sigma==sigmaprime); % Identity matrix
            MPO_Even{sigma,site,sigmaprime} = MatrixStorageRE{sigma,2*(mod(site,2)==1)+1*(mod(site,2)==0),sigmaprime}; % Right edge Hamiltonian
        end
    end   
end
clear('MatrixStorageB'); clear('MatrixStorageLE'); clear('MatrixStorageRE');
end