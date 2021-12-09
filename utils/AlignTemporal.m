function phase = AlignTemporal( phase, wrap_index, weights )
%AlignTemporal Temporally align displacement phases after unwrapping
% phase        the unwrapped phase images (dim: Nx * Ny * ..)
% wrap_index   the dimension of the phase offsets 
% (weights)    optional: Weights to guide alignment process (dim: Nx * Ny)
%              e.g. use averaged magnitude over all dynamics & phases
%
% This function "aligns time" by finding the most suitable phase wrap combination that temporally
% unwraps most of the data.
%
% Author: Christian Guenthner, 21/10/2019 - ETH Zurich, guenthner@biomed.ee.ethz.ch 
    
    % reorder phase
    order = [wrap_index setdiff(1:ndims(phase),wrap_index)];
    phase = permute(phase,order);     
    sz = size(phase);
    
	% unwrap all phase images in time:
	cor = round((phase - unwrap(phase,[],1))/pi/2);
	cor = reshape(cor,sz(1),sz(2)*sz(3),[]);
	
	% if weights are not supplied, use ones instead
	if nargin<3 || isempty(weights)
		weights = ones([1 sz(1:end)]);
		weights = reshape(weights,sz(1),sz(2)*sz(3),[]);
	else
		weights = permute(weights,order);
		weights = reshape(weights,sz(1),sz(2)*sz(3),[]);
	end
	assert(all(size(weights)==size(cor)),'Weights and phase must be of equal size')
	%
	unicor = unique(cor(~isnan(cor(:))));
	
	cor = permute(cor,[2 1 3]);
	weights = permute(weights,[2 1 3]);
	
	M = weightedmode(cor,weights,unicor);
	
	cor = reshape(M,sz(1),1,size(cor,3));
	cor = cor * 2 *pi;

	% bring phase into same size as correction matrix
    phase = reshape(phase,sz(1),sz(2)*sz(3),[]);    
	% align phase
    phase = phase - cor;
    
    % reshape  & backpermute
    phase = reshape(phase,sz);
    phase = ipermute(phase,order);
end

function M=weightedmode(S,W,cat)
    S=S(:,:);
    W=W(:,:);
    assert(all(size(W)==size(S)),'Data and weights must be equally large');
    assert(all(isreal(S(:))),'Data must be real');
    assert(all(isreal(W(:))),'Weights must be real');

    if nargin<3 || isempty(cat)
        cat=unique(S(:));
    else
        cat=cat(:);
    end
    
    sz = size(S);
    
    csum=zeros(length(cat),sz(2));
    
    for n=1:length(cat)
        csum(n,:) = sum(W .* (S==cat(n)),1,'omitnan');
    end
    
    [~,I] = max(csum,[],1);    
    M = reshape(cat(I),1,[]);
end
