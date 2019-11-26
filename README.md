# MATLAB-Dataspace-to-Figure-Units
Convert from dataspace to figure units to make it easier to add annotations pointing to data in a MATLAB figure window.

The annotation function, which allows you to programmatically add a wide range of annotations to your figure, requires coordinates to be specified in normalized figure units. I have found that I almost always want to specify my annotations in data space (i.e., based on the values of data displayed in an axes).
This utility function converts coordinates in data space into normalized figure coordinates, for input to annotation. Some annotations require you to specify (x,y) pairs, while others require a 4 element position vector. This function supports both syntaxes.

Here's a simple example:

```matlab
% Create some data
t = 0:.1:4*pi;
s = sin(t);

% Add an annotation requiring (x,y) coordinate vectors
plot(t,s);ylim([-1.2 1.2])
xa = [1.6 2]*pi; % X-Coordinates in data space
ya = [0 0]; % Y-Coordinates in data space
[xaf,yaf] = ds2nfu(xa,ya); % Convert to normalized figure units
annotation('arrow',xaf,yaf) % Add annotation
```

Note: I believe annotation was introduced in MATLAB 7.

[![View Data space to figure units conversion on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion)
