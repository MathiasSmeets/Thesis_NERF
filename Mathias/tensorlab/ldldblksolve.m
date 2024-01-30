function X = ldldblksolve(L,D,b)
%LDLKRONSOLVE short discription
%   long description
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/04/08   NV      Initial version
    
    if size(L,1) == size(L,2)
        % If not given in compact format, convert to compact format
        I = find(L(:,1),2);
        if numel(I) > 1, I = I(2)-1; 
        else I = size(L,1); end
        R = size(L,1)/I;

        i = repmat((1:R*I).',1,R);
        j = bsxfun(@plus, repmat((1:I).',R,1), (0:R-1)*I);
        L = L((j - 1)*R*I + i);
    else 
        R = size(L,2);
        I = size(L,1)/R;
    end 

    if nargin < 3
        % Extract D from L
        b = D;
        D = L(ones(I,1) * ((0:R-1)*I*R) + reshape(1:R*I,I,R));
    end 
    
    % If D not in compact format, convert
    if numel(D) == (R*I)^2, D = diag(D); end
    
    L = reshape(L, I, R, R);
    D = reshape(D, I, R);
    X = reshape(b, I, R);
    
    tol = max(abs(D(:))) * eps(class(D)) * numel(D) * 10;

    % Solve x = L\b
    for r = 2:R
        X(:,r) = (X(:,r) - sum(squeeze(L(:,r,1:r-1)) .* X(:,1:r-1), 2)); 
    end 
    % Solve x = D\x (using pseudoinverse)
    X = X ./ D;
    X(abs(D) < tol) = 0; 
    % Solve x = L'\x
    for r = R-1:-1:1
        X(:,r) = X(:,r) - sum(conj(L(:,r+1:R,r)) .* X(:,r+1:R), 2);
    end 
    % Fix output format depending on b
    X = reshape(X, size(b));
end
