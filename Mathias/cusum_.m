function [iUpper, iLower, upperSum, lowerSum] = cusum(x, varargin)
%CUSUM  Detect small changes in mean via cumulative sums
% [IUPPER, ILOWER] = CUSUM(X) returns the first index of the upper and
% lower cumulative sums of vector X that have moved beyond 5 standard
% deviations above and below the target mean, respectively.  IUPPER and
% ILOWER will be empty if all indices are within the tolerance.  The
% minimum detectable mean shift to detect is set to one standard deviation.
% The target mean and standard deviations are estimated from the first 25
% samples in the input signal.
%
% [IUPPER, ILOWER] = CUSUM(X, CLIMIT) sets the control limits that the
% upper and lower cumulative sums are allowed to drift from the mean.
% Specify the control limit, CLIMIT, in units of standard deviations from
% the mean. The default value of CLIMIT is 5.
%
% [IUPPER, ILOWER] = CUSUM(X, CLIMIT, MSHIFT) sets the minimum mean shift
% to detect.  Specify the mean shift, MSHIFT, in units of standard
% deviations away from the mean.  The default value of MSHIFT is 1.
%
% [IUPPER, ILOWER] = CUSUM(X, CLIMIT, MSHIFT, TMEAN) specifies the target
% mean of the overall signal from which to make the baseline measurement.
% If unspecified, TMEAN is computed as the mean of the first 25 samples of
% X.
%
% [IUPPER, ILOWER] = CUSUM(X, CLIMIT, MSHIFT, TMEAN, TDEV) specifies the
% target standard deviation from which to compute the upper and lower
% control limits.  If unspecified, TDEV is computed as the first 25 samples
% of X.
%
% [IUPPER, ILOWER] = CUSUM(...,'all') will return all indices of the upper
% and lower cumulative sums that are beyond the control limit.
%
% [IUPPER, ILOWER, UPPERSUM, LOWERSUM] = CUSUM(...) additionally returns
% the upper and lower cumulative sums.
%
% CUSUM(...) without output arguments plots the upper and lower
% cumulative sums normalized to one standard deviation above and below the
% target mean, respectively.
%
%   % Example:
%   %    Detect change in mean of a signal that has a mean shift in
%   %    its data at the 50th sample point the signal.
%   load('meanshift','x')
%   plot(x)
%   title('Original Signal')
%   figure
%   cusum(x)
%
%   See also MEAN, FINDCHANGEPTS.

%   Copyright 2015-2020 The MathWorks, Inc.
%#codegen

narginchk(1,6);

if nargin > 1
    inputArgs = cell(1,length(varargin));
    [inputArgs{:}] = convertStringsToChars(varargin{:});
else
    inputArgs = varargin;
end

% validate input vector
validateattributes(x,{'numeric'},{'vector','real','finite'},'cusum','X',1);

if iscolumn(x)
    rowFlag = false;
else
    rowFlag = true;
end

xCol = x(:);

% get presence of 'all' input option and check for any unrecognized string
[allLimits,inputArgs1] = getAllOption(inputArgs{:});

% fetch and validate numeric arguments
[climit,mshift,tmean,tdev] = getNumericArguments(xCol,inputArgs1{:});

% compute the upper and lower cumulative sums
[uppersum,lowersum] = cusums(xCol,tmean,mshift,tdev);

% fetch the indices where the sums exceed the borders
[iupper,ilower] = violations(uppersum,lowersum,climit,tdev,allLimits);

if nargout == 0
    coder.internal.assert(coder.target('MATLAB'),'signal:cusum:PlottingNotSupported');
    % plot control chart
    plotchart(uppersum,lowersum,tmean,tdev,climit,iupper,ilower);
else
    % copy variables to output
    if rowFlag
        iUpper = iupper.';
        iLower = ilower.';
        upperSum = uppersum.';
        lowerSum = lowersum.';
    else
        iUpper = iupper;
        iLower = ilower;
        upperSum = uppersum;
        lowerSum = lowersum;
    end
end

end
% -------------------------------------------------------------------------
function [allOption,argsOut] = getAllOption(varargin)

allOption1 = false;

coder.unroll();
for i = 1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
        coder.internal.assert(coder.internal.isConst(arg), ...
            'signal:cusum:FlagAsConst');
        coder.internal.assert(strncmpi(arg,'all',length(arg)), ...
            'signal:cusum:UnrecognizedOption',arg);
        allOption1 = true;
        idx = i;
    end
end

if allOption1
    allOption = true;
    argsOut = cell(1,length(varargin)-1);
    [argsOut{:}] = varargin{[1:(idx-1),(idx+1):end]};
else
    allOption = false;
    argsOut = varargin;
end

end
% -------------------------------------------------------------------------
function [climit,mshift,tmean,tdev] = getNumericArguments(x,varargin)

nvarargin = numel(varargin);

% input data type
xDataType = class(x);
% samples to compute mean and standard deviation
xSamples = x(1:min(25,length(x)));

if nvarargin < 4
    tdev = max(eps(xDataType),std(xSamples));
else
    tdev1 = varargin{4};
    validateattributes(tdev1,{'numeric'},{'scalar','real','positive','finite'},'cusum','TDEV',5);
    tdev = cast(tdev1(1),xDataType);
end

if nvarargin < 3
    tmean = mean(xSamples);
else
    tmean1 = varargin{3};
    validateattributes(tmean1,{'numeric'},{'scalar','real','finite'},'cusum','TMEAN',4);
    tmean = cast(tmean1(1),xDataType);
end

if nvarargin < 2
    mshift = 1;
else
    mshift1 = varargin{2};
    validateattributes(mshift1,{'numeric'},{'scalar','real','positive','finite'},'cusum','MSHIFT',3);
    mshift = double(mshift1(1));
end

if nvarargin < 1
    climit = 5;
else
    climit1 = varargin{1};
    validateattributes(climit1,{'numeric'},{'scalar','real','positive','finite'},'cusum','CLIMIT',2);
    climit = double(climit1(1));
end

end
% -------------------------------------------------------------------------
function [uppersum,lowersum] = cusums(x,tmean,mshift,tdev)

uppersum = coder.nullcopy(zeros(size(x),'like',x));
lowersum = coder.nullcopy(zeros(size(x),'like',x));

uppersum(1) = 0;
lowersum(1) = 0;

for i = 2:length(x)
    uppersum(i) = max(0, uppersum(i-1) + x(i) - tmean - mshift*tdev/2);
    lowersum(i) = min(0, lowersum(i-1) + x(i) - tmean + mshift*tdev/2);
end

end
% -------------------------------------------------------------------------
function [iupper,ilower] = violations(uppersum,lowersum,climit,tdev,allLimits)

if allLimits
    iupper = find(uppersum >  climit*tdev);
    ilower = find(lowersum < -climit*tdev);
else
    iupper = find(uppersum >  climit*tdev, 1, 'first');
    ilower = find(lowersum < -climit*tdev, 1, 'first');
end

end
% -------------------------------------------------------------------------
function  plotchart(uppersum,lowersum,tmean,tdev,climit,iupper,ilower)
newplot
n = length(uppersum);
hLines = plot(1:numel(uppersum),uppersum(:)/tdev,1:numel(lowersum),lowersum(:)/tdev);
line([1 n],climit*[1 1],'LineStyle',':','Color',hLines(1).Color);
line([1 n],-climit*[1 1],'LineStyle',':','Color',hLines(2).Color);

if ~isempty(iupper)
    line(iupper,uppersum(iupper)/tdev,'Color','r','Marker','o','LineStyle','none');
end

if ~isempty(ilower)
    line(ilower,lowersum(ilower)/tdev,'Color','r','Marker','o','LineStyle','none');
end

hAxes = hLines(1).Parent;
hAxes.YLim = [min(-climit-1,hAxes.YLim(1)) max(climit+1,hAxes.YLim(2))];
ylabel(getString(message('signal:cusum:StandardErrors')))
xlabel(getString(message('signal:cusum:Samples')))
title(getString(message('signal:cusum:CUSUMControlChart', ...
    sprintf('\n\\mu_{target} = %f',tmean), ...
    sprintf('\\sigma_{target} = %f',tdev))));
end