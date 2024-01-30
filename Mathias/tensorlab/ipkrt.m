function M = ipkrt(U, varargin)
%IPKRT Inner product of transposed Khatri-Rao matrix with itself.
%   IPKRT(A,B) computes the inner product krt(A,B)'*krt(A,B) in which krt is
%   the transposed or row-wise Khatri-Rao product. Blocking is used to reduce
%   memory overhead.
%      
%   IPKRT(A,B,C,...) and IPKRT({A B C ...}) compute the inner product M'*M of a
%   string of transposed Khatri-Rao products M = A x B x C x ..., where x
%   denotes the transposed Khatri-Rao product.
%  
%   IPKRT(A,B,'BlockSize',b) controls the maximal block size, i.e., the number
%   of rows in each matrix A, B, used to construct M to limit memory usage.
%   Default = 1e9/8/prod(size_U), in which size_U = cellfun('size',
%   {A,B,...}, 2).
%
%   IPKRT(A,B,'Subs',subs) and IPKRT(U,'Subs',subs) with subs a cell of
%   length(U) containing row sampling indices for U, i.e., M = U{1}(subs{1},:) x
%   U{2}(subs{2},:) x ..., where x denotes the transposed/row-wise Khatri-Rao
%   product.
%
%   See also krt, inprod, ipkrt2.
    
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
    
    if nargin >= 2
        % Split factors and options
        factidx = find(cellfun(@ischar, varargin),1);
        if isempty(factidx)
            factors = varargin;
            varargin = {};
        else 
            factors = varargin(1:factidx-1);
            varargin = varargin(factidx:end);
        end 
    else 
        factors = [];
    end 
    
    % Collect factors
    if ~iscell(U), U = [{U} factors];
    else U = [U factors]; end
    size_U = cellfun('size', U, 2);
    
    % Parse options
    p = inputParser;
    p.addOptional('Subs', []);
    p.addOptional('BlockSize', ceil(1e9/8/prod(size_U)));
    p.parse(varargin{:});
    p.KeepUnmatched = true;
    p.PartialMatching = true;
    options = p.Results;
    
    % Check selection indices
    if ~isempty(options.Subs)
        if ~iscell(options.Subs) || numel(options.Subs) ~= numel(U)
            error('ipkrt:subs',['options.Subs should be a cell with a length ' ...
                                'equal to the number of input matrices.']);
        end
        options.Subs = cellfun(@(s) s(:), options.Subs, 'UniformOutput', ...
                               false);
        I = length(options.Subs{1});
        if any(cellfun('length', options.Subs) ~= I)
            error('ipkrt:subs',['numel(options.Subs{n}) should be equal for ' ...
                                'all n.']);
        end 
    else 
        I = size(U{end}, 1);
        if any(cellfun('size',U,1) ~= I)
            error('ipkrt:U','Input matrices should have the same number of rows.');
        end
    end

    % Compute first block
    i = min(I, options.BlockSize);
    if isempty(options.Subs)
        Ui = cellfun(@(u) u(1:i,:), U, 'UniformOutput', false);
    else
        Ui = cellfun(@(u,s) u(s(1:i),:), U, options.Subs, 'UniformOutput', ...
                     false);
    end    
    tmp = krt(Ui);
    M = tmp'*tmp;
    
    % Compute remaining blocks
    while i < I
        idx = i+1:min(i+options.BlockSize,I);
        if isempty(options.Subs)
            Ui = cellfun(@(u) u(idx,:), U, 'UniformOutput', false);
        else 
            Ui = cellfun(@(u,s) u(s(idx),:), U, options.Subs, 'UniformOutput', ...
                         false);
        end 
        tmp = krt(Ui);
        M = M + tmp'*tmp;
        i = i + options.BlockSize;
    end 
end
