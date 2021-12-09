function data = DecodePhase( data, encoding_scheme, dim_dynamic )
%DECODEPHASE Decode encoding directions to actual displacement field
% data             the unwrapped and aligned phase images (dim: Nx * Ny * ..)
% encoding_scheme  either PMSR or Hadamard
% dim_dynamic      dimension of the different encodings 
%
% This function takes unwrapped & temporally aligned phases and decodes the encoded displacement
% fields by calculating the pseudo-inverse of the gradient encoding matrix. 
%
% For a detailed description of encoding and decoding in MRE, see:
% 
%    Guenthner, C., & Kozerke, S. (2018). Encoding and readout strategies in magnetic resonance elastography. 
%    NMR in Biomedicine, 31(10), e3919. https://doi.org/10.1002/nbm.3919
%
% Author: Christian Guenthner, 21/10/2019 - ETH Zurich, guenthner@biomed.ee.ethz.ch 
    
    % get encoding matrix for classical PMSR encoding
    enc_matrix = [1 0 0 0; 0 1 0 0; 0 0 1 0; 1 1 1 1].';
	
	% if hadamard, change it accordingly
    if strcmpi(encoding_scheme,'Hadamard')
        enc_matrix = [1 -1 1 -1;1 1 -1 -1; 1 -1 -1 1; 1 1 1 1].';
    elseif any(strcmpi(encoding_scheme,{'PMSR','PMS'}))
        % do nothing. enc_matrix is already set.
    else
        error('Encoding scheme `%s` is not supported.',encoding_scheme);
    end
    
    % decoding operation is given by pseudo-inverse
    dec_matrix = pinv(enc_matrix);
    
    % permute dynamic index to front
    indxs = [dim_dynamic setdiff(1:length(size(data)),dim_dynamic)];
    data = permute(data,indxs);
    
    % get size of data & reshape
    sz = size(data);    
    data = reshape(data,size(dec_matrix,2),[]);
	
	%decoding operation
    data = dec_matrix * data;

	% remove phase offset dimension
    data = data(1:3,:);
    sz(1) = 3;
    
    % reshape data
    data = reshape(data,sz);
    %permute data back
    data = ipermute(data,indxs);   
end
