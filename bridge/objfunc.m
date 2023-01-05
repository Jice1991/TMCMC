function  y=objfunc(x,y)
npar = y.n; % N0. parameters
freqtrue = y.freq; 
modeltrue = y.mode;
s1 = y.s1;
s2 = y.s2;
Ns=20;nn=5;mm=5;
%% Load Sensitivity matrix
alpha_act = x';
LoadStructure
n_modes = 5; % Number of measured modes
modeIndex = 1 : n_modes; % Indexes of these measured modes

%% Assemble structure matrices
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;

%% Simulate "experimental data"
[psiExpAll,lambdaExp] = eigs(K_act,M0,n_modes,'sm');
[lambdaExp,dummyInd] = sort((diag(sqrt(lambdaExp)/2/pi)),'ascend') ;
lambdaExp = lambdaExp(modeIndex);
psiExpAll = psiExpAll(:,dummyInd(modeIndex));
psi_m = psiExpAll(measDOFs,:); % 21 DOFs

% Normalize the mode shape vectors by maximum entry
for i = 1: n_modes
    psi_m(:,i) = psi_m(:,i) / norm(psi_m(:,i));
end
freq = lambdaExp(1:n_modes);
phi = psi_m(:,1:n_modes);

freq1=repmat(freq,1,Ns);
D1=(freqtrue-freq1).^2; % Original Code

phi=remodel(phi);
for i=1:mm
    for j=1:Ns
        a=modeltrue(:,i,j);
        b=phi(:,i);
        D2(i,j)= (a-b)'*(a-b)./s2(i); % Original Code

    end
end
 
 y=-Ns*log(2*pi)-Ns/2*log(sum(s1))-Ns/2*log(sum(s2))-0.5*(sum(sum(D1')./s1) +  sum(sum(D2')));  % Original Code



