function waveField = gradwrapFFT(smoothedPhase, inplaneResolution, parameters)

% get the dimensions
n1 = size(smoothedPhase,1);% number of rows
n2 = size(smoothedPhase,2);% number of columns
nSlice = size(smoothedPhase,3);% number of slices
nTimestep = size(smoothedPhase,4);% number of time steps
nComponent = size(smoothedPhase,5);% number of components
nFrequency = size(smoothedPhase,6);% number of frequencies

% loop for unwrapping and fft
gradient1 = zeros(n1, n2, nSlice, nTimestep);
gradient2 = zeros(n1, n2, nSlice, nTimestep);
waveField = zeros(n1, n2, nSlice, 2, nComponent, nFrequency);
for iFrequency = 1 : nFrequency
    for iComponent = 1 : nComponent
        
        for iTimestep = 1 : nTimestep
            currentData = smoothedPhase(:,:,:,iTimestep,iComponent,iFrequency);
            
            if nSlice > 1
                [temp2, temp1] = gradient(exp(1i*currentData), inplaneResolution(2), inplaneResolution(1), 1);
            else % nSlice == 1
                [temp2, temp1] = gradient(exp(1i*currentData), inplaneResolution(2), inplaneResolution(1));
            end
            
            % in-plane derivative components
            temp = exp(-1i*currentData);
            gradient1(:,:,:,iTimestep) = imag(temp1 .* temp);
            gradient2(:,:,:,iTimestep) = imag(temp2 .* temp);
        end
        
        % fourier transformation and selection of harmonic
        fourier1 = fft(gradient1,[],4);
        waveField(:,:,:,1,iComponent,iFrequency) = fourier1(:,:,:,1+parameters.numberOfHarmonics);
        fourier2 = fft(gradient2,[],4);
        waveField(:,:,:,2,iComponent,iFrequency) = fourier2(:,:,:,1+parameters.numberOfHarmonics);
        
    end
end

end