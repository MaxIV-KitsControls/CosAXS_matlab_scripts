function [Iout, varargout] = normalize( q, I, normrange, varargin )
% Oskar Berntsson, 2021
% Iout = normalize( q, I, normrange )
% normalizes the data in I column-wise, using only the specified q-range.

if numel(normrange) == size(I,2)
    % in this case normrange is values provided to be normalized with.
    % Could be I0, or could be something else.
    Iout = I ./ repmat(normrange,numel(q),1);
else    
    normMask = normrange(1)<=q & normrange(2)>q;
    Iout = I ./ repmat(mean(I(normMask,:), 1), length(q), 1);
end


end

