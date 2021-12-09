% in-plane Laplacian unwrap
%
% [unwrapWaveField] = gunwrap(waveField)
%
% input:
% waveField                 6D temporal resoluted wave field
%                               - 1st and 2nd dimensions correspond to the in-plane matrix
%                                 size(number of rows and columns)
%                               - 4th dimension are the number of timesteps
% 
% output:
% unwrapWaveField          6D unwraped wave field

% edit by Florian Dittmann, Charité-Universitätsmedizin Berlin, Berlin, Germany
% idea from Ingolf Sack and Florian Dittmann, Charité-Universitätsmedizin Berlin, Berlin, Germany
% Date: 2015/04/01
% last change: 2018/10/09 (by Heiko Tzschätzsch)
function [unwrapWaveField] = gunwrap(waveField)

% size of waveField
si = size(waveField);
si = [si 1 1 1 1]; % length(si)>=6

% Dirichlet condition: displacement at the edge = 0;
E = ones(si(1)*si(2),1);
laplaceMatrix = spdiags([E E -4*E E E], [-si(1) -1:1 si(1)], si(1)*si(2) , si(1)*si(2) );
laplaceMatrix(si(1):si(1):end, si(1)+1:si(1):end) = 0;% Dirichlet for bottom row
laplaceMatrix(si(1)+1:si(1):end, si(1):si(1):end) = 0;% Dirichlet for top row

unwrapWaveField = zeros(si);
si(3) = 1;
for iSlice = 1 : size(waveField,3)
    % rewrite exp(1i*phase) of  all single 2D phase images as column vectors in a matrix
    normMRSignal = reshape(exp(1i*waveField(:,:,iSlice,:,:,:)),[si(1)*si(2) prod(si(3:end))]);
    
    % compute 2D laplacian of phase images
    laplacianField = imag((laplaceMatrix*normMRSignal).*conj(normMRSignal));
    
    % solve integration
    currentComplexWaveField = laplaceMatrix\laplacianField;
    
    unwrapWaveField(:,:,iSlice,:,:,:) = reshape(currentComplexWaveField, si);
end

end