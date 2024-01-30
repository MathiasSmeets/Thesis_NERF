%# Template for CPD GUI code generator
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
% CPD computations 
% 
% File generated on ${timestamp} using GUI_CPD.

%% Load data 
$IF{fromfile}
tmp       = load('${file}','${var}');
T         = tmp.${var};       
$END
N         = getorder(${tensor});
size_tens = getsize(${tensor});

%% Set data options
R = ${R};
$IF{usesymmetry}
% Symmetry definition
variablesinmode = ${varinmode};
$END
$IF{usenonneg}
% Nonnegativity definition
nonnegmodes     = ${nonnegmodes};
$END


%% Set algorithm options
options = struct;
options.Display          = 10;    	      
options.TolFun           = ${TolFun};    
options.TolX             = ${TolX};
options.MaxIter          = ${MaxIter};
options.Algorithm        = ${algorithm}; 
$IF{showcurves}          
options.ShowCurves       = true;
$END
$IF{usenonneg}
options.NonnegativeModes = nonnegmodes;
$END
$IF{usesymmetry}
options.VariablesInMode  = variablesinmode;
$END
$IF{autoinit}
options.Initialization   = ${initalgorithm};
$END
$IF{randinit}
options.InitializationOptions = struct('Real', ${dist});
$END


%% Call the cpd function
[${U}, output] = cpd(${tensor}, ${Uinit}, options);

%% Compute metrics/properties of result
% Results contains the various metrics
results = analyze_cpd(${tensor}, ${U}, 'ShowPlots', false);

%% Generate plots
$IF{showcurves}
figure('Name', 'Convergence curves');
plot_convergence(output);
$END

figure('Name', 'Angles between rank-1 terms');
plot_angles(${U});
title('Angles between rank-1 terms');

figure('Name', 'Error per slice');
plot_errslice(${tensor}, ${U});
title('Error per slice');

figure('Name', 'Computed CPD versus original data');
visualize(${U}, 'original', ${tensor});

figure('Name', 'Norms of rank-1 terms'); 
stem(results.norms, 'filled')
title('Norms of rank-1 terms');

plot_rank1terms(${U});

