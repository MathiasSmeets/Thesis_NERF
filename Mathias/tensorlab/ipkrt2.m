function M = ipkrt2(U1, U2, varargin)
%IPKRT2 Inner product of two transposed Khatri-Rao matrices.
%   IPKRT2(U1,U2) computes the inner product krt(U1{:})'*krt(U2{:}) in which krt
%   is the transposed or row-wise Khatri-Rao product. Blocking is used to reduce
%   memory overhead.
%    
%   See also krt, inprod, ipkrt.
    
%   Authors: Nico Vervliet       (Nico.Vervliet@kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
%   Version History:
%   - 2017/04/04   NV      Documentation and tests
%   - 2016/08/18   NV      Initial version
% 
%   References: 
%   [1] Vervliet, N., Debals, O., De Lathauwer, L., "Canonical polyadic
%       decomposition of incomplete tensors with linearly constrained factors",
%       Technical Report 16-172, ESAT-STADIUS, KU Leuven, Leuven, Belgium, 2017.
    
    size_U1 = cellfun('size', U1, 2);
    size_U2 = cellfun('size', U2, 2);
    
    % Parse options
    p = inputParser;
    p.addOptional('BlockSize', ceil(1e9/8/(prod(size_U1)+prod(size_U2))));
    p.parse(varargin{:});
    p.KeepUnmatched = true;
    p.PartialMatching = true;
    options = p.Results;
    
    I = size(U1{end}, 1);
    if any(cellfun('size',U1,1) ~= I) || any(cellfun('size',U2,1) ~= I)
        error('ipkrt2:U','Input matrices U1{n}, U2{n} should have the same number of rows.');
    end
    
    % First block
    i = min(I, options.BlockSize);
    U1i = cellfun(@(u) u(1:i,:), U1, 'UniformOutput', false);
    U2i = cellfun(@(u) u(1:i,:), U2, 'UniformOutput', false);
    tmp1 = krt(U1i);
    tmp2 = krt(U2i);
    M = tmp1'*tmp2;
    % Remaining blocks
    while i < I
        idx = i+1:min(i+options.BlockSize,I);
        U1i = cellfun(@(u) u(idx,:), U1, 'UniformOutput', false);
        U2i = cellfun(@(u) u(idx,:), U2, 'UniformOutput', false);
        tmp1 = krt(U1i);
        tmp2 = krt(U2i);
        M = M + tmp1'*tmp2;
        i = i + options.BlockSize;
    end 
end
