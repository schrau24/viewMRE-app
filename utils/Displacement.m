function displacement = Displacement( data, fft_index)
%Displacement Apply fourier transformation on data to obtain displacements
% data          unwrapped, aligned and decoded phase data
% fft_index     wave phase offset dimension
%
% This function applies the Fourier transformation along the wave phase dimension to 
% obtain the complex displacement vector field.
%
% Author: Christian Guenthner, 21/10/2019 - ETH Zurich, guenthner@biomed.ee.ethz.ch 

    % permute dynamic index to front
    indxs = [fft_index setdiff(1:length(size(data)),fft_index)];
    data = permute(data,indxs);
    
    % get size of data & reshape
    sz = size(data);    
    datap = reshape(data,sz(1),[]);

    % number of phase offsets
    N=size(datap,1);    
    
	% let's put together an inverse FFT operator:
	[X,Y]=meshgrid((0:(N-1)),(0:(N-1)));
	iDFT = exp(1i.*X.*Y/N.*2.*pi);
	DFT = pinv(iDFT);

    % and apply:
	data = DFT*datap*2i;
    
    
    % Get Displacementmap (first frequency component)
    displacement = data(2,:); %[radians]
    
	% Restructure Data, Reshape and iPermute
    sz(1) = 1;
	displacement = reshape(displacement,sz);
    displacement = ipermute(displacement,indxs);
end