function varargout = ds2nfu(varargin)
% DS2NFU  Convert data space units into normalized figure units. 
%
% [Xf, Yf] = DS2NFU(X, Y) converts X,Y coordinates from
% data space to normalized figure units, using the current axes.  This is
% useful as input for ANNOTATION.  
%
% POSf = DS2NFU(POS) converts 4-element position vector, POS from
% data space to normalized figure units, using the current axes.  The
% position vector has the form [Xo Yo Width Height], as defined here:
%
%      web(['jar:file:D:/Applications/MATLAB/R2006a/help/techdoc/' ...
%           'help.jar!/creating_plots/axes_pr4.html'], '-helpbrowser')
%
% [Xf, Yf] = DS2NFU(HAX, X, Y) converts X,Y coordinates from
% data space to normalized figure units, on specified axes HAX.  
%
% POSf = DS2NFU(HAX, POS) converts 4-element position vector, POS from
% data space to normalized figure units, using the current axes. 
%
% Ex.
%       % Create some data
% 		t = 0:.1:4*pi;
% 		s = sin(t);
%
%       % Add an annotation requiring (x,y) coordinate vectors
% 		plot(t,s);ylim([-1.2 1.2])
% 		xa = [1.6 2]*pi;
% 		ya = [0 0];
% 		[xaf,yaf] = ds2nfu(xa,ya);
% 		annotation('arrow',xaf,yaf)
%
%       % Add an annotation requiring a position vector
% 		pose = [4*pi/2 .9 pi .2];
% 		posef = ds2nfu(pose);
% 		annotation('ellipse',posef)
%
%       % Add annotations on a figure with multiple axes
% 		figure;
% 		hAx1 = subplot(211);
% 		plot(t,s);ylim([-1.2 1.2])
% 		hAx2 = subplot(212);
% 		plot(t,-s);ylim([-1.2 1.2])
% 		[xaf,yaf] = ds2nfu(hAx1,xa,ya);
% 		annotation('arrow',xaf,yaf)
% 		pose = [4*pi/2 -1.1 pi .2];
% 		posef = ds2nfu(hAx2,pose);
% 		annotation('ellipse',posef)

% Michelle Hirsch
% mhirsch@mathworks.com
% Copyright 2006-2014 The MathWorks, Inc

%% Process inputs
narginchk(1,3)

% Determine if axes handle is specified
if isscalar(varargin{1}) && ishandle(varargin{1}) && strcmp(get(varargin{1},'type'),'axes')	
	hAx = varargin{1};
	varargin = varargin(2:end);
else
	hAx = gca;
end;

errmsg = ['Invalid input.  Coordinates must be specified as 1 four-element \n' ...
	'position vector or 2 equal length (x,y) vectors.'];

% Proceed with remaining inputs
if isscalar(varargin)	% Must be 4 elt POS vector
	pos = varargin{1};
	if length(pos) ~=4
		error(errmsg);
    end
else
	[x,y] = deal(varargin{:});
	if length(x) ~= length(y)
		error(errmsg)
	end
end

	
%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');
if hAx.DataAspectRatioMode == "manual" || hAx.PlotBoxAspectRatioMode == "manual"
    axpos = plotboxpos(hAx);
else
    axpos = get(hAx,'Position');
end

% When DataAspectRatioMode and PlotBoxAspectRatioMode are both manual,
% Position likely won't accurately reflect actual position of visible axis

axlim = axis(hAx);

if strcmp(hAx.XScale,"log")
    axlim(1:2) = log10(axlim(1:2)); 
	x = log10(x);
end
if strcmp(hAx.YScale,"log")
    axlim(3:4) = log10(axlim(3:4)); 
	y = log10(y);
end


axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));

%% Transform data
if exist('x','var')
	varargout{1} = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
	varargout{2} = (y-axlim(3))*axpos(4)/axheight + axpos(2);
else
	pos(1) = (pos(1)-axlim(1))/axwidth*axpos(3) + axpos(1);
	pos(2) = (pos(2)-axlim(3))/axheight*axpos(4) + axpos(2);
	pos(3) = pos(3)*axpos(3)/axwidth;
	pos(4) = pos(4)*axpos(4)/axheight;
	varargout{1} = pos;
end


%% Restore axes units
set(hAx,'Units',axun)
end

function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
axisPos  = getpixelposition(h);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    xlim = get(h, 'XLim');
    ylim = get(h, 'YLim');
    
    % Deal with axis limits auto-set via Inf/-Inf use
    
    if any(isinf([xlim ylim]))
        hc = get(h, 'Children');
        hc(~arrayfun( @(h) isprop(h, 'XData' ) & isprop(h, 'YData' ), hc)) = [];
        xdata = get(hc, 'XData');
        if iscell(xdata)
            xdata = cellfun(@(x) x(:), xdata, 'uni', 0);
            xdata = cat(1, xdata{:});
        end
        ydata = get(hc, 'YData');
        if iscell(ydata)
            ydata = cellfun(@(x) x(:), ydata, 'uni', 0);
            ydata = cat(1, ydata{:});
        end
        isplotted = ~isinf(xdata) & ~isnan(xdata) & ...
                    ~isinf(ydata) & ~isnan(ydata);
        xdata = xdata(isplotted);
        ydata = ydata(isplotted);
        if isempty(xdata)
            xdata = [0 1];
        end
        if isempty(ydata)
            ydata = [0 1];
        end
        if isinf(xlim(1))
            xlim(1) = min(xdata);
        end
        if isinf(xlim(2))
            xlim(2) = max(xdata);
        end
        if isinf(ylim(1))
            ylim(1) = min(ydata);
        end
        if isinf(ylim(2))
            ylim(2) = max(ydata);
        end
    end

    dx = diff(xlim);
    dy = diff(ylim);
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis

hparent = get(h, 'parent');
hfig = ancestor(hparent, 'figure'); % in case in panel or similar
currax = get(hfig, 'currentaxes');

temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', hparent);
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);

set(hfig, 'currentaxes', currax);
end




