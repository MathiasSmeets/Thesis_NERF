%# Template for multilinear compression GUI code generator
%# 
%# Lines starting with $# are not exported. 
%# Variables are denoted by ${varname}.
%# The following conditional statements are allowed 
%#
%#     $IF{condition} true code $ELSE false code $END
%#
%#     $IF{condition}
%#     true code  (multiline possible)
%#     $ELSE 
%#     false code (multiline possible)
%#     $END
%# 
%# The condition refers to a field in the conditionals struct, which evaluates to true or
%# false. If tests cannot be nested. Else statements can be omitted. 
%#
% Multilinear compression 
% 
% File generated on ${timestamp} using GUI_MLSVD.

%% Load data 
$IF{fromfile}
tmp       = load('${file}','${var}');
T         = tmp.${var};       
$END
N         = getorder(${tensor});
size_tens = getsize(${tensor});

%% Compute compression
size_core = ${sizecore};
[${U},${S},sv]  = ${algorithm}(${tensor}, size_core);

$IF{refinement}
% Refine result
[${U},${S}] = ${refinement}(${tensor}, ${U}, ${S});
$END

%% Visualize residual
figure;
visualize({${U},${S}}, 'original', ${tensor});

%% Compute errors and compression rate
err_abs = froblmlrares(${tensor}, ${U}, ${S});
err_rel = err_abs/frob(${tensor});
fit_rel = 1 - err_rel;

% compression ratio
ratio   = numel(${tensor}) / (sum(cellfun(@numel,${U}))+numel(${S}));

%% Find entries with maximal error
% Absolute error
residual = lmlrares(${tensor}, ${U}, ${S});
[max_err_abs, ind] = max(abs(residual(:)));
sub = cell(1,N);
[sub{:}] = ind2sub(size_tens, ind);
ind_err_abs = horzcat(sub{:});
fprintf('Maximal absolute error of %6.3e at %s.\r\n', max_err_abs, mat2str(ind_err_abs));

% Relative error
[max_err_abs, ind] = max(abs(residual(:)./${tensor}(:)));
[sub{:}] = ind2sub(size_tens, ind);
ind_err_rel = horzcat(sub{:});
fprintf('Maximal relative error of %6.3e at %s.\r\n', max_err_abs, mat2str(ind_err_abs));

%% Plot errors per slice
plot_errslice(${tensor}, ${U}, ${S}, 'absolute')
