function [absG, phi, strain] = laplaceInversion(shearWaveField, inplaneResolution, frequency, parameters)

% get the dimensions
n1 = size(shearWaveField,1);% number of rows
n2 = size(shearWaveField,2);% number of columns
nSlice = size(shearWaveField,3);% number of slices
nGrad = size(shearWaveField,4);% number of gradients
nComponent = size(shearWaveField,5);% number of components
nFrequency = size(shearWaveField,6);% number of frequencies

% loop for MDEV-inversion
numeratorPhi = zeros(n1,n2,nSlice);
denominatorPhi = zeros(n1,n2,nSlice);
numeratorAbsG = zeros(n1,n2,nSlice);
denominatorAbsG = zeros(n1,n2,nSlice);
strain = zeros(n1,n2,nSlice);
for iFrequency = 1 : nFrequency
    for iComponent = 1 : nComponent
        for iGrad = 1 : nGrad
            
            U = shearWaveField(:, :, :, iGrad, iComponent, iFrequency);
            
            if nSlice > 1
                [gradient2, gradient1] = gradient(U,inplaneResolution(2), inplaneResolution(1), 1);
                [~, gradient11] = gradient(gradient1,inplaneResolution(2), inplaneResolution(1), 1);
                [gradient22, ~] = gradient(gradient2,inplaneResolution(2), inplaneResolution(1), 1);
            else % nSlice == 1
                [gradient2, gradient1] = gradient(U,inplaneResolution(2), inplaneResolution(1));
                [~, gradient11] = gradient(gradient1,inplaneResolution(2), inplaneResolution(1));
                [gradient22, ~] = gradient(gradient2,inplaneResolution(2), inplaneResolution(1));
            end
            
            % laplace of U
            LaplaceU = gradient11 + gradient22;
            
            % calculation of numerator and denominator for MDEV-inversion
            numeratorPhi = numeratorPhi + real(LaplaceU).*real(U) + imag(LaplaceU).*imag(U);
            denominatorPhi = denominatorPhi + abs(LaplaceU).*abs(U);
            
            numeratorAbsG = numeratorAbsG + parameters.density*(2*pi*frequency(iFrequency)).^2.*abs(U);
            denominatorAbsG = denominatorAbsG + abs(LaplaceU);
            
            % sum the strain
            strain = strain + abs(U);
            
        end
    end
end

% avoid division by zero
denominatorPhi(denominatorPhi == 0) = eps;
denominatorAbsG(denominatorAbsG == 0) = eps;

% inversion
phi = acos(-numeratorPhi./denominatorPhi); %[rad]
if isfield(parameters,'noiseCorrection') && parameters.noiseCorrection
    phi = phi - 0.2;
end
absG = numeratorAbsG./denominatorAbsG; %[Pa]

end