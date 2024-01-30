%CPD_JACOBIAN_INCOMPLETE CPD Jacobian of incomplete tensor.
%  CPD_JACOBIAN_INCOMPLETE(J, U, SUB) performs the in-place computation of the
%  (transposed) Jacobian J for a CPD with factor matrices U{1}, ..., U{N} of an
%  incomplete Nth-order tensor T, i.e., each entry J(i,k) is the derivative of
%  the kth known entry of CPDGEN(U) - T w.r.t. the ith variable. The
%  (multilinear) indices or subscripts of the known entries are given in the
%  cell sub such that the kth known entry has multilinear index (sub{1}(k), ...,
%  sub{N}(k)). The values of J are changed in-place, i.e., the non-zero pattern
%  is not changed. U{n} should all be either real or complex, stored as a
%  double. U{n} should contain the transposed factor matrices. Sub should be an
%  int32 cell and the order the known entries should match the pattern of J.
%
%  Inputs:
%  - J         Sparse matrix corresponding to transposed Jacobian,
%              initialized with the correct pattern. If U{n} are real (complex), J
%              should be real (complex). Should be a double.
%  - U         Transposed factor matrices, i.e., the rank equals the number
%              of rows each factor matrix U{n}. Should be a double. 
%  - sub       Cell of multilinear indices. Should be int32. Length(sub) can
%              be smaller than length(U) if size(U{n},2) = 1 for n >
%              length(sub). 
%  Outputs:
%  - J         Updated transposed Jacobian.
%
%  Usage:
%  
%      cpd_jacobian_incomplete(J, U, T.sub); % No output argument!
%
%  Compilation: This file should be compiled first using mex:
%
%      mex -largeArrayDims cpd_jacobian_incomplete.c

% Author(s):  Nico Vervliet       (Nico.Vervliet@kuleuven.be)
%             Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2017/04/04   NV      Complex data and order(U) > order(T) added.
% - 2017/03/06   NV      Initial version.
