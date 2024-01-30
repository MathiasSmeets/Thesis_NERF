function X = krt(U,varargin)
%KRT Transposed or row-wise Khatri-Rao product.
%   KRT(A,B) returns the transposed or row-wise Khatri-Rao product of two
%   matrices A and B, of dimensions I-by-J and I-by-K respectively. The result
%   is an I-by-J*K matrix formed by the matching row-wise Kronecker products,
%   i.e., the kth row of the transposed Khatri-Rao product is defined as
%   kron(A(k,:),B(k,:)).
%
%   KR(A,B,C,...) and KR({A B C ...}) compute a string of transposed Khatri-Rao
%   products A x B x C x ..., where x denotes the transposed Khatri-Rao product.
%
%   See also kron, kr, ipkrt, ipkrt2.

%   Authors: Nico Vervliet       (Nico.Vervliet@kuleuven.be)
%            Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)

if ~iscell(U), U = [{U} varargin];
else, U = [U varargin]; end

[I,K] = size(U{end});
if any(cellfun('size',U,1) ~= I)
    error('krt:U','Input matrices should have the same number of rows.');
end

X = U{end};
for n = length(U)-1:-1:1
    J = size(U{n},2);
    A = reshape(U{n},[I 1 J]);
    X = reshape(bsxfun(@times,A,X),[I J*K]);
    K = J*K;
end
X = reshape(X,[I, size(X,2)]);
