%maxInjection
%
%  [mscores, binds] = maxInjection( scoresAB )
%
%  Find maximum score mapping from elements A to B.  Uses Matlab's linear
%  programming function in the optimization toolbox.
%  scoresAB: [Na x Nb] score for mapping A to B
%          Na <= Nb
%  binds: [1 x Na] gives column index (of B) for each A
%                  
%  Daniel Morris, Feb 2018
%
%  See also: maxSegDice
%

%  Copyright (c) 2018, Daniel Morris
%  Licensed under the Simplified BSD License
function [mscores, binds] = maxInjection( scoresAB )

[nr,nc] = size(scoresAB);
if nr > nc, error('Na must be <= Nb for an injection'); end

%want 1 selection per row (gt)
%want at most 1 selection per column (segment)
colcon = arrayfun(@(x) sparse(ones(1,nr), nr*(x-1)+(1:nr), ones(1,nr), 1, nr*nc), 1:nc, 'uniformoutput', false );
rowcon = arrayfun(@(x) sparse(ones(1,nc), (x-1)+1:nr:nr*nc, ones(1,nc), 1, nr*nc), 1:nr, 'uniformoutput', false );

A = vertcat( vertcat(colcon{:}), vertcat(rowcon{:}) );

%Do linear programming to maximize dice:
x = linprog( -scoresAB(:), A, ones(nr+nc,1), [], [], zeros(nr*nc,1), ones(nr*nc,1), optimoptions('linprog','Display','None') );

x = reshape(x, nr, nc );

if sum(x>0) > nr
    warning('Not an injection');
end

[~, binds] = max(x, [], 2);
binds = binds';

mscores = scoresAB( sub2ind([nr nc],1:nr, binds) );

end

