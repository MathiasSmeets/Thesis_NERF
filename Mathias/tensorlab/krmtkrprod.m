function M = krmtkrprod(T,Ul,Ur,m,n)
%KRMTKRPROD Khatri-Rao times matricized tensor times Khatri-Rao product.
%   KRMTKRPROD(T,UL,UR,m,n) matricizes the tensor and multiplies from the left
%   with the conjugated tranpose of the Khatri-Rao product of all matrices in UL
%   except the mth, and from the right with the Khatri-Rao product of all
%   matrices in UR expect the nth. At most two permutations of smaller
%   matrices are needed. In inefficient matlab code, the product is as
%   follows:
%   
%       tmp = kr(UL([end:-1:m+1 m-1:-1:1]))' * ...
%                 tens2mat(T,[1:m-1 m+1:length(UL)]);
%       tmp = reshape(tmp, [size(tmp,1) size(T,m) cellfun('size',UR,1)]);
%       tmp = tens2mat(tmp, [2 1 n+2])*kr(UR([end:-1:n+1 n-1:-1:1]));
%       M   = reshape(tmp, size(T,m)*RL, size(T,n+3)*RR)
%
%   RL and RR are the number of columns of the factor matrices in UL and UR,
%   respectively. Note the conjugated transpose on UL.
%
%   Note: inputs are unchecked for performance reasons.
%    
%   See also mtkrprod.
    
%   Authors: Nico Vervliet       (Nico.Vervliet@kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
%   Version History:
%   - 2017/04/04   NV      Documentation and tests
%   - 2016/08/10   NV      Initial version
% 
%   References: 
%   [1] Vervliet, N., Debals, O., De Lathauwer, L., "Canonical polyadic
%       decomposition of incomplete tensors with linearly constrained factors",
%       Technical Report 16-172, ESAT-STADIUS, KU Leuven, Leuven, Belgium, 2017.
    
    Nl = length(Ul);
    Nr = length(Ur);
    U  = [Ul, Ur];
    n  = n  + Nl;
    N  = Nr + Nl;
   
    Rl = size(Ul{1}, 2);
    Rr = size(Ur{1}, 2);
    sz = getsize(T);
    currsz = sz;
    
    % First contractions on free indices using matrix products
    if m > 1
        M = reshape(T, sz(1), prod(currsz(2:end)));
        u = Ul{1};
        M = u'*M;
        currsz(1) = Rl;
        if n < N
            M = reshape(M, prod(currsz(1:end-1)), sz(N));
            u = Ur{end};
            M = M*u;
            currsz(end) = Rr;
        end 
    elseif n < N
        M = reshape(T, [], sz(N))*Ur{end};
        currsz(end) = Rr;
    else 
        M = T;
    end 
    
    % Left-to-right contractions
    if m >= 3, M = reshape(M, currsz(1), prod(currsz(2:end))); end 
    for k = 2:m-1
        tmp = M;
        M = zeros(Rl, prod(currsz(3:end)));
        for r = 1:Rl
            tmp2 = U{k}(:,r)'*reshape(tmp(r,:),currsz(2),[]);
            M(r,:) = tmp2(:);
        end 
        currsz(2) = [];
    end 
    % Right-to-left contractions
    if n <= N-2, M = reshape(M, prod(currsz(1:end-1)), currsz(end)); end 
    for k = N-1:-1:n+1
        tmp = M;
        M = zeros(size(M,1)/sz(k),Rr);
        for r = 1:Rr
            tmp2 = reshape(tmp(:,r),[],sz(k))*U{k}(:,r);
            M(:,r) = tmp2(:);
        end 
        currsz(end-1) = [];
    end

    % Compute permutation of data
    if m == 1, permleft = [2:Nl, 1];
    else permleft = [1, 3:Nl-m+2, 2]; end 
    if n == N, permright = [n-Nl, 1:Nr-1];
    else permright = [n-Nl, 1:n-Nl-1, n-Nl+1]; end 
    perm = [permleft, length(permleft)+permright];

    % Permute to middle
    M = reshape(M, currsz);
    M = permute(M, perm);
    currsz = currsz(perm);
    
    % Use matrix multiplications if possible
    if m == 1
        M = reshape(M, currsz(1), prod(currsz(2:end)));
        M = Ul{2}'*M;
        currsz(1) = Rl;
    end 
    if n == N
        M = reshape(M, prod(currsz(1:end-1)), currsz(end));
        M = M*Ur{end-1};        
        currsz(end) = Rr;
    end 
        
    % Use left-to-right contractions for remaining matrix in UL
    M = reshape(M, Rl, prod(currsz(2:end)));
    for k = m+1+(m==1):Nl
        tmp = M;
        M = zeros(Rl, prod(currsz(3:end)));
        for r = 1:Rl
            tmp2 = U{k}(:,r)'*reshape(tmp(r,:),currsz(2),[]);
            M(r,:) = tmp2(:);
        end 
        currsz(2) = [];
    end 
    % Use right-to-left contractions for remaining matrix in UR
    M = reshape(M, prod(currsz(1:end-1)),Rr);
    for k = n-1-(n==N):-1:Nl+1
        tmp = M;
        M = zeros(prod(currsz(1:end-2)),Rr);
        for r = 1:Rr
            tmp2 = reshape(tmp(:,r),[],currsz(end-1))*U{k}(:,r);
            M(:,r) = tmp2(:);
        end 
        currsz(end-1) = [];
    end
    
    % Permute end result
    M = permute(reshape(M, Rl, sz(m), sz(n), Rr), [2 1 3 4]);
    M = reshape(M, sz(m)*Rl, sz(n)*Rr);
end
