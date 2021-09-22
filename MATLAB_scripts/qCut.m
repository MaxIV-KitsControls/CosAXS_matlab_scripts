function [ I, q ] = qCut( I, q, qrange )
% Oskar Berntsson, 2021
% [ I, q ] = qCut( I, q, qrange )
% removes the rows if I that do not fall in the specified q range. This
% changes the dimension of I, and so a new q vector is returned.

I = I(q >= min(qrange) & q <= max(qrange), :);
q = q(q >= min(qrange) & q <= max(qrange));


end

