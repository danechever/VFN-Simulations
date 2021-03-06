function [ n, m ] = zern_noll2index( noll_index )
%zern_noll2index Summary of this function goes here
%   Detailed explanation goes here

% Convert from noll index to Zernike indices
n = 0;
j1 = noll_index-1;
while(j1 > n)
    n = n + 1;
    j1 = j1 - n;
end

m = (-1)^noll_index  * (mod(n,2) + 2 * floor((j1+mod(n+1,2)) / 2.0 ));

end

