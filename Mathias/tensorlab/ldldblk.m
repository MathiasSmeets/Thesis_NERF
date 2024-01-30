function [A,D] = ldldblk(A,compressed)
%LDLDBLK Compute LDL factorization of diagonal block matrix.
%   [L,D] = LDLDBLK(A) computes the LDL factorization of A, i.e., A = L*D*L', in the case A is an
%   I*J x I*J Hermitian matrix in which each of the J*J blocks is a I x I diagonal matrix, e.g., A =
%   kron(W,eye(I)) with W an J x J Hermitian matrix. Both L and D are given in a compressed format,
%   i.e., only the diagonals are returned in as an I*J x J matrix and an I*J vector,
%   respectively. A can be given in the same format as well.
%
%   L = LDLDBLK(A) stores the values of D on the diagonal of L.
%
%   [L,D] = LDLDBLK(A,false) returns the full rather than the compressed format.

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/07/25   NV      Complex data and uncompressed format.
% - 2018/04/08   NV      Initial version
    
    if nargin < 2, compressed = true; end
    
    if size(A,1) == size(A,2)
        % If not given in compact format, convert to compact format
        I = find(A(:,1),2);
        if numel(I) > 1, I = I(2)-1; 
        else I = size(A,1); end
        R = size(A,1)/I;

        i = repmat((1:R*I).',1,R);
        j = bsxfun(@plus, repmat((1:I).',R,1), (0:R-1)*I);
        A = A((j-1)*R*I + i);
    else 
        R = size(A,2);
        I = size(A,1)/R;
    end 
    
    A = reshape(A,I,R,R); % L is stored in A.
    D = zeros(I,R);
    for r = 1:R
        if r == 1
            D(:,r) = A(:,1,1);
        else 
            v = zeros(I, 1, r-1);
            for s = 1:r-1
                v(:,1,s) = bsxfun(@times, conj(A(:,r,s)), A(:,s,s));
            end  
            D(:,r) = A(:,r,r) -  sum(A(:,r,1:r-1) .* v ,3);
        end 
        A(:,1:r-1,r) = 0;
        zidx = abs(D(:,r)) > eps(class(D));
        
        if sum(zidx) > 0
            if r == 1
                A(:,r+1:R,r) = A(:,r+1:R,r)./D(:,r);
            else
                A(:,r+1:R,r) = (A(:,r+1:R,r) - sum(A(:,r+1:R,1:r-1).*v,3))./D(:,r);
            end
        end
        A(~zidx,r+1:R,r) = 0;
        A(:,r,r) = D(:,r);
    end    
    if nargout > 1
        for r = 1:R
            A(:,r,r) = 1;
        end 
    end
    A = reshape(A,I*R,R);
    D = reshape(D,I*R,1);
    
    if ~compressed
        i = repmat((1:R*I).',1,R);
        j = bsxfun(@plus, repmat((1:I).',R,1), (0:R-1)*I);
        A = full(sparse(i,j,A));
        D = diag(D);
    end 
end
