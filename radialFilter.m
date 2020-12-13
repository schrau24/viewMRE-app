
function shearWaveField = radialFilter(waveField, inplaneResolution, parameters)

% get the dimensions
n1 = size(waveField,1);% number of rows
n2 = size(waveField,2);% number of columns
nSlice = size(waveField,3);% number of slices
nGrad = size(waveField,4);% number of gradients
nComponent = size(waveField,5);% number of components
nFrequency = size(waveField,6);% number of frequencies

% calculate the wavenumber for k-space filtering
k1 = -( (0:n1-1)-fix(n1/2) ) / (n1*inplaneResolution(1));%[1/m] wavenumber in 1st direction
k2 = ( (0:n2-1)-fix(n2/2) ) / (n2*inplaneResolution(2));%[1/m] wavenumber in 2nd direction
absK = hypot( ones(n1,1)*k2, k1'*ones(1,n2) );%[1/m] transform into polar coordinates

filter = 1 ./ (1 + (absK/parameters.lowpassThreshold).^(2*parameters.lowpassOrder));
filter = ifftshift(filter);

% loop for filtering
shearWaveField = zeros(size(waveField));
for iFrequency = 1 : nFrequency
    for iComponent = 1 : nComponent
        for iGrad = 1 : nGrad
            for iSlice = 1 : nSlice
                
                currentWaveField = waveField(:, :, iSlice, iGrad, iComponent, iFrequency);
                
                % filtering in k-space
                data = fftn(currentWaveField);
                filteredData = data .* filter;
                currentShearWaveField = ifftn(filteredData);
                
                shearWaveField(:, :, iSlice, iGrad, iComponent, iFrequency) = currentShearWaveField;
                
            end
        end
    end
end

end
