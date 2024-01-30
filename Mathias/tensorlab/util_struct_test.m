function util_struct_test(z,f,fz,relerr,verbose)
%UTIL_STRUCT_TEST Utility for tests of structures
%   UTIL_STRUCT_TEST(z,f,fz) with z the structure variables, f(z,task) the
%   structure function, and fz the intended function evaluation for z,
%   compares the function evaluation and the left and right
%   Jacobian-vector products (using DERIV).
%
%   UTIL_STRUCT_TEST(...,relerr) uses relative errors relerr(1) for the
%   comparison of fz and f(z,[]), and uses relerr(2) and relerr(3) for the
%   left and right Jacobian-vector products respectively.
%    (default [1e-7 1e-3 1e-3])

%   Author: Otto Debals         (Otto.Debals@esat.kuleuven.be)
%           Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%           Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
%   Version History:
%   - 2016/09/01    NV      Added size test
%   - 2014/03/05    OD      Initial version
%   - 2014/03/25    OD      Cell variables possible with serialize

   
if nargin < 4 || isempty(relerr), relerr = [1e-7 1e-3 1e-3]; end
if nargin < 5, verbose = false; end
if verbose, display = @(s) disp(s); 
else display = @(s) 0; end

display('Computing function value');
[x,task] = f(z,[]);
assertElementsAlmostEqual(fz,x,'relative',relerr(1));

display('Computing Jacobian');
J = deriv(@(z) f(z,[]),z,[],'Jacobian');

display('Computing right Jacobian-vector product');
try 
    task.l = []; task.r = [];
    task.r = creaternd(z);
    Jr = J*serialize(task.r);
    Jr2 = f(z,task);
    assertElementsAlmostEqual(Jr(:),Jr2(:),'relative',relerr(2));
    assertElementsAlmostEqual(size(Jr2), size(fz));
catch e
    if ~isempty(strfind(e.identifier, 'nonanalytic'))
        display('Not analytic. Skipping.')
    else 
        rethrow(e)
    end 
end 


[~,task] = f(z,[]);
task.l = []; task.r = [];
task.l = creaternd(fz);

display('Computing left Jacobian-vector product');
Jl = J'*task.l(:);
Jl2 = f(z,task);
Jl2s = serialize(Jl2);
assertElementsAlmostEqual(Jl(:),Jl2s(:),'relative',relerr(3));
if iscell(Jl2)
    sz1 = cellfun(@size, Jl2, 'UniformOutput', false);
    sz2 = cellfun(@size, z, 'UniformOutput', false);
    assertEqual(size(sz1), size(sz2));
    for k = 1:length(sz1)
        assertElementsAlmostEqual(sz1{k}, sz2{k});
    end 
else 
    assertElementsAlmostEqual(size(Jl2), size(z));
end 

end

function x = creaternd(z)
% Creates an x having the same structure as z but with randomized content
% Supports structures, cells, logicals, characters and numerics (real and
% complex)

if isstruct(z)
    x = z;
    flds = fieldnames(z);
    for i = 1:numel(flds)
        x.(flds{i}) = creaternd(z.(flds{i}));
    end
elseif iscell(z)
    x = cellfun(@(y) creaternd(y),z,'UniformOutput',false);
elseif islogical(z)
    x = rand(size(z))>0.5;
elseif ischar(z)
    symbols = ['a':'z' 'A':'Z' '0':'9'];
    nums = randi(numel(symbols),size(z));
    x = symbols(nums);
elseif isnumeric(z)
    if isreal(z), x = randn(size(z));
    else x = randn(size(z))+1i*randn(size(z));
    end
else
    error('util_struct_test:creaternd:type','Type is not supported!')
end

end

function z = serialize(z)
% Serializes the data in the cell z into a vector

    if iscell(z)
        for i = find(cellfun(@iscell,z(:).'))
            z{i} = serialize(z{i});
        end
        s = cellfun(@numel,z(:)); o = [0; cumsum(s)];
        c = z; z = zeros(o(end),1);
        for i = 1:length(s), z(o(i)+(1:s(i))) = c{i}(:); end
    else
        z = z(:);
    end
end