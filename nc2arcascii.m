function nc2arcascii(file, var, dimOrder, outputbase, varargin)
%NC2ARCASCII Converts netCDF file to Arc ascii file
%
% nc2arcascii(file, var, dimOrder, outputbase)
% nc2arcascii(file, var, dimOrder, outputbase, depthIndices)
% nc2arcascii(file, var, dimOrder, outputbase, depthIndices, timeformat)
% nc2arcascii(file, var, dimOrder, outputbase, depthIndices, timeformat,...
%             cellsize)
%
% This function converts a netCDF file into the Arc ascii grid format
% requested by the UBC team.  It has been developed with GFDL-style netCDF
% files in mind; it may or may not work on all other files.  This function
% takes the specified variable, breaks it into separate lat x lon grids
% (one for each time and specified depth), regrids to a square-cell grid if
% necessary, recenters the grid on the Atlantic Ocean, and prints each grid
% to a new file.
%
% Input variables:
%
%   file:           name of .nc file
%
%   var:            name of variable to be added to new files
%
%   dimOrder:       cell array telling the dimensions of the variable.
%                   Each cell must hold one of the following strings:
%                   'time', 'depth', 'lat', 'lon'.
%
%   outputbase:     base name for output file
%
%   depthIndices:   indices of depth slices you want to extract.  If not
%                   included and the variable in question is 4D, all depths
%                   will be extracted (if variable is 3D, this is ignored).
%                   An index of -1 indicates that you want to extract the
%                   deepest data points, i.e. the data along the ocean
%                   bottom.
%
%   timeformat:     format used for units of time variable.  Include if not
%                   a typical datestr format.  If not included, Matlab will
%                   interpret the datestring as explained in datestr
%                   documentation.
%
%   cellsize:       size of cells in ouput grid, in degrees.  If not
%                   included, default is to use the finest resolution found
%                   in the netCDF file.
%
% Example:
%
%   I want to convert advection data at the surface and the bottom of the 
%   ocean to Arc ascii files.  The data is currently held in the netCDF
%   file example.nc as a variable named u.  For this example, I chose a
%   more complicated scenario that requires all 7 inputs; depending on the
%   file, you may need anywhere from 4 to 7 inputs.
%
%   > ncdump -h example.nc
%
%     netcdf example {
%     dimensions:
%         time = UNLIMITED ; // (1200 currently)
%         nv = 2 ;
%         xu_ocean = 360 ;
%         yu_ocean = 200 ;
%         zt_ocean = 50 ;
%         zt_edges_ocean = 51 ;
%     variables:
%         double time(time) ;
%             time:long_name = "time" ;
%             time:units = "days since 2001-01-01 00:00:00" ;
%             time:cartesian_axis = "T" ;
%             time:calendar_type = "NOLEAP" ;
%             time:calendar = "NOLEAP" ;
%             time:bounds = "time_bounds" ;
%         float nv(nv) ;
%             nv:long_name = "vertex number" ;
%             nv:units = "none" ;
%             nv:cartesian_axis = "N" ;
%         float xu_ocean(xu_ocean) ;
%             xu_ocean:long_name = "ucell longitude" ;
%             xu_ocean:units = "degrees_E" ;
%             xu_ocean:cartesian_axis = "X" ;
%         float yu_ocean(yu_ocean) ;
%             yu_ocean:long_name = "ucell latitude" ;
%             yu_ocean:units = "degrees_N" ;
%             yu_ocean:cartesian_axis = "Y" ;
%         float zt_ocean(zt_ocean) ;
%             zt_ocean:long_name = "tcell depth" ;
%             zt_ocean:units = "meters" ;
%             zt_ocean:cartesian_axis = "Z" ;
%             zt_ocean:positive = "down" ;
%             zt_ocean:edges = "zt_edges_ocean" ;
%         float zt_edges_ocean(zt_edges_ocean) ;
%             zt_edges_ocean:long_name = "tcell depth edges" ;
%             zt_edges_ocean:units = "meters" ;
%             zt_edges_ocean:cartesian_axis = "Z" ;
%             zt_edges_ocean:positive = "down" ;
%         float u(time, zt_ocean, yu_ocean, xu_ocean) ;
%             u:long_name = "zonal current" ;
%             u:units = "m/sec" ;
%             u:valid_range = -10.f, 10.f ;
%             u:missing_value = -10.f ;
%             u:cell_methods = "time: mean" ;
%             u:time_avg_info = "average_T1,average_T2,average_DT" ;
%         double average_T1(time) ;
%             average_T1:long_name = "Start time for average period" ;
%             average_T1:units = "days since 2001-01-01 00:00:00" ;
%         double average_T2(time) ;
%             average_T2:long_name = "End time for average period" ;
%             average_T2:units = "days since 2001-01-01 00:00:00" ;
%         double average_DT(time) ;
%             average_DT:long_name = "Length of average period" ;
%             average_DT:units = "days" ;
%         double time_bounds(time, nv) ;
%             time_bounds:long_name = "time axis boundaries" ;
%             time_bounds:units = "days" ;
% 
%     // global attributes:
%             :filename = "ocean.220101-230012.u.nc" ;
%             :title = "CM2.1U-H2_SresA1B_X1" ;
%     }
%
%
%   The first two inputs, for file and variable name, are straightforward:
%   'example.nc' and 'u', repectively.  The description shows that u is a
%   variable u(time, zt_ocean, yu_ocean, xu_ocean), which I know is time,
%   depth, latitude, and longitude, respectively, so the third input to the
%   function will be {'time', 'depth', 'lat', 'lon'}.  The fourth input is
%   the basename for output files; I choose 'example'.  Because I only want
%   the surface (depth #1) and along the bottom, I add a fifth input
%   specifying this: [1 -1].  In this file, the date format for the time
%   variable units is not one Matlab natively recognizes (if you're not
%   sure whether a format is acceptable, you can try to run nc2arcascii
%   with the defaults; it will  error out with an explanatory message if
%   the format is unacceptable).  So the sixth input is 'yyyy-mm-dd
%   HH:MM:SS'.  Finally, I want the output grids to have a 1 degree
%   resolution, so I include 1 as the final input.
%
%   My final command would therefore be:
%
%   nc2arcascii('example.nc', 'u', {'time','depth','lat','lon'}, ...
%               'example', [1 -1], 'yyyy-mm-dd HH:MM:SS', 1);
%
%   This will produce 2400 output files, one for each of the 1200 times
%   stored in the file and for each of the two depths I requested.  They
%   will be named example_yyyymmdd_depthX.asc, where yyyymmdd is the date
%   assciated with each time slice, as determined by the values of the time
%   variable in the file, and X is the depth slice index (-1 results in
%   'bottom' rather than depthX).

% Copyright 2007 Kelly Kearney

%-------------------------
% Determine dimensions
%-------------------------

fprintf('\nnc2arcascii:\n%s\n\n', file);

% Find variable of interest

FileInfo = nc_info(file);
varIndex = find(strcmp(var, {FileInfo.DataSet.Name}));

if isempty(varIndex)
    error('Specified variable not found in file');
end

% Dimension names

[tf, dimLoc] = ismember({'time', 'depth', 'lat', 'lon'}, dimOrder);

latName = FileInfo.DataSet(varIndex).Dimension{dimLoc(3)};
lonName = FileInfo.DataSet(varIndex).Dimension{dimLoc(4)};
timeName = FileInfo.DataSet(varIndex).Dimension{dimLoc(1)};

hasDepth = tf(2);
if hasDepth
    depthName = FileInfo.DataSet(varIndex).Dimension{dimLoc(2)};
end

% Check if time boundary data variable is present

possibleNames = {sprintf('%s_bounds', timeName), ...
                 sprintf('%s_bnds', timeName), ...
                 sprintf('%s_bnd', timeName)};

[tf, tbndIndex] = ismember(possibleNames, {FileInfo.DataSet.Name});
if any(tf)
    hastbnd = true;
    tbndIndex = tbndIndex(tf);
    tbndName = FileInfo.DataSet(tbndIndex).Name;
else
    hastbnd = false;
end

%-------------------------
% Parse optional 
% timeformat and cellsize
% inputs
%------------------------

if nargin > 5 && ~isempty(varargin{2})
   timeformat = varargin{2};
else
    timeformat = '';
end

if nargin > 6
    cellSize = varargin{3};
else
    cellSize = NaN;
end


%-------------------------
% Read lat, lon, time, and 
% depth data
%-------------------------

fprintf('Reading lat, lon, time, and depth data...\n');

latdata  = nc_varget(file, latName);
londata  = nc_varget(file, lonName);
if hastbnd
    timedata = nc_varget(file, tbndName);
else
    timedata = nc_varget(file, timeName);
end
    
nlat = length(latdata);
nlon = length(londata);
ntime = size(timedata,1);

if hasDepth
    
    depthData = nc_varget(file, depthName);
    
    % If not specified, extract all depths
    
    if nargin < 5
        depthIndices = 1:length(depthdata);
    else
        depthIndices = varargin{1};
    end
   
    ndepth = length(depthIndices); 
    
else
    ndepth = 1;
end

%-------------------------
% Calculations to see if
% and how grids will 
% need to be adjusted
%-------------------------

fprintf('Regridding if necessary...\n');

% Check if original grid is uniform with square cells

minLatCell = min(diff(latdata));
minLonCell = min(diff(londata));

latIsUniform = all(diff(latdata) == minLatCell);
lonIsUniform = all(diff(londata) == minLonCell);

isSquare = latIsUniform && lonIsUniform && (minLatCell == minLonCell);

% If no cell size is provided, use finest resolution found in file

if isnan(cellSize)
    cellSize = min([minLatCell minLonCell]);
end

needsRegrid = ~(isSquare & cellSize == minLatCell);

% Regrid lat and lon if necessary

if ~needsRegrid
    londataNew = londata;
    latdataNew = latdata;
else
    [latdataGrid, londataGrid] = meshgrat(latdata,londata);
    [temp, refvec] = geoloc2grid(latdataGrid, londataGrid, zeros(nlat,nlon), cellSize);
    [latgrat, longrat] = meshgrat(temp, refvec); % data irrelevant, only for size
    latdataNew = latgrat(:,1);
    londataNew = longrat(1,:);
end
nlatNew = length(latdataNew);
nlonNew = length(londataNew);

% Index order to center grid on Atlantic (order lon from -180 to 180)
       
londataNew(londataNew > 180)  = londataNew(londataNew > 180) - 360;
londataNew(londataNew < -180) = londataNew(londataNew < -180) + 360;
[londataNew, isortlon] = sort(londataNew);

% Index order in order to put northermost values at top of grid

[latdataNew, isortlat] = sort(latdataNew, 1, 'descend');

%-------------------------
% Parse time and depth 
% info.  This will be used 
% to name the files later, 
% but taking care of it
% here so problems can be 
% caught before the time-
% consuming stuff
%-------------------------

fprintf('Checking time units...\n');

% Get time strings associated with each time (I assume that the time
% variable has an attribute for units, which has the format "months[or
% days] since datestr")

timeIndex = find(strcmp(timeName, {FileInfo.DataSet.Name})); 
[tf, unitLoc] = ismember('units', {FileInfo.DataSet(timeIndex).Attribute.Name});
timeUnit = FileInfo.DataSet(timeIndex).Attribute(unitLoc).Value;

[timeStep, temp1, temp2] = strread(timeUnit, '%s since %s%s', 1);
timeStep = timeStep{1};
datestring = [temp1{1} ' ' temp2{1}];
if isempty(timeformat)
    try
        startDate = datevec(datestring);
    catch
        error('Unrecognized date format in time variable units\n\n%s\n\nPlease give date format as input variable', datestring);
    end
else
    startDate = datevec(datestring, timeformat);
end

timeStr = arrayfun(@(x) futuredate(startDate, x, timeStep), timedata, 'uni', 0);


if size(timeStr,2) == 2 % time bounds
    timeStr = cellfun(@(x,y) sprintf('%s_to_%s', x,y), timeStr(:,1), timeStr(:,2), 'uni', 0);
end

% timeVec = repmat(startDate, length(timedata), 1);
% switch timeStep
%     case 'days'
%         timeVec(:,3) = timeVec(:,3) + timedata;
%     case 'months'
%         timeVec(:,2) = timeVec(:,2) + timedata;
%     otherwise
%         error('This file uses an unknown time unit for step size\n\n%s', timeStep);
% end
% 
% timeStr = cellstr(datestr(datenum(timeVec), 'yyyymmdd'));

% Get depth string associated with each depth.  Because varying t-cell
% depth may be confusing if used in a filename, only the depth index rather
% than the actual depth is used. 

if hasDepth
    depthStr = cell(ndepth,1);
    for idepth = 1:ndepth
        if depthIndices(idepth) == -1
            depthStr{idepth} = 'bottom';
        else
            depthStr{idepth} = sprintf('depth%d', depthIndices(idepth));
        end
    end 
end

%-------------------------
% Read and process 
% variable of interest at 
% each time and write to 
% ascii file
%-------------------------

fprintf('Extracting main variable and writing to file...\n');

for itime = 1:ntime
    
    fprintf('  Time %d of %d\n', itime, ntime);
    
    % For crashed runs only, skips iterations that have already been
    % completed
    
    if nargin > 7 && strcmp(varargin{4}, 'crashmode')
        for idepth = 1:ndepth
            if hasDepth
                testfile{idepth} = sprintf('%s_%s_%s_%s.asc', outputbase, var, timeStr{itime}, depthStr{idepth});
            else
                testfile{idepth} = sprintf('%s_%s_%s.asc', outputbase, var, timeStr{itime});
            end
        end
        alreadyran = cellfun(@(x) exist(x, 'file'), testfile);
        if all(alreadyran)
            continue
        end
    end
    
    % Read one time slice from file
   
    if hasDepth
        
        % If file holds 4D array, data is 3D depth x lat x lon
        
        startInd = zeros(1,4);
        startInd(dimLoc(1)) = itime-1;
        stepInd = repmat(-1,1,4);
        stepInd(dimLoc(1)) = 1;
        
        data = nc_varget(file, var, startInd, stepInd);
        
    else
        
        % If file holds 3D array, data is 2D lat x lon
        
        startInd = zeros(1,3);
        startInd(dimLoc(1)) = itime-1;
        stepInd = repmat(-1,1,3);
        stepInd(dimLoc(1)) = 1;
        
        data = nc_varget(file, var, startInd, stepInd);
     
    end
    
   
    for idepth = 1:ndepth
        
        fprintf('    Depth %d of %d\n', idepth, ndepth);
        
        % Extract one depth slice from data array if it is 3D, or just save
        % if it is 2D
        % TODO Right now assumes lat x lon are the trailing dimensions.
        % Add some permutations so it can handle any dimesion order
        
        if hasDepth
            if depthIndices(idepth) == -1
                
                % -1 indicates to get deepest points, i.e. ocean bottom 
                 
                hasdata = ~isnan(data);
                ibottom = squeeze(sum(hasdata));
                [ilon, ilat] = meshgrid(1:nlon,1:nlat);

                vardata = NaN(nlat,nlon);

                island = ibottom == 0;
                ind = sub2ind(size(data), ibottom(~island), ilat(~island), ilon(~island));
                vardata(~island) = data(ind);

            else
                
                % Other indices represent a flat slice from the 3D array
                 
                vardata = squeeze(data(depthIndices(idepth),:,:));
            end
      
        else
            vardata = data;
        end

        
        % Regrid if necessary
        
        if needsRegrid
            [vardata, refvec] = geoloc2grid(latdataGrid, londataGrid, vardata, cellSize);
        end      
        
        vardata = vardata(isortlat,isortlon);
        
        % Replace missing values with -9999
        
        [tf, missingValLoc] = ismember('missing_value', {FileInfo.DataSet(varIndex).Attribute.Name});
        if tf
            missingVal = FileInfo.DataSet(varIndex).Attribute(missingValLoc).Value;
            isMissing = vardata == missingVal;
            vardata(isMissing) = -9999;
        end
        
        % Create file name for new file
        
        if hasDepth
            filename = sprintf('%s_%s_%s_%s.asc', outputbase, var, timeStr{itime}, depthStr{idepth});
        else
            filename = sprintf('%s_%s_%s.asc', outputbase, var, timeStr{itime});
        end
        
        % Write header
        
        fid = fopen(filename, 'wt');
        
        fprintf(fid, 'ncols %d\n', nlonNew);
        fprintf(fid, 'nrows %d\n', nlatNew);
        fprintf(fid, 'xllcorner %f\n', min(londataNew));
        fprintf(fid, 'yllcorner %f\n', min(latdataNew));
        fprintf(fid, 'cellsize %f\n', cellSize);
        fprintf(fid, 'nodata_value -9999\n');
        
        fclose(fid);
        
        % Write data
        
        dlmwrite(filename, vardata, '-append', 'delimiter', ' ');
        
    end
end

fprintf('Done\n');


%------------------------------
% Subfunction: Project
% date string using calendar
% with no leap year
%------------------------------

function str = futuredate(start, timepassed, unit)

switch unit
    
    case 'months'
        
        yearsPassed = floor(timepassed/12);
        monthsPassed = mod(timepassed,12);
        
        newDate = start + [yearsPassed monthsPassed 0 0 0 0];
        str = datestr(newDate, 'yyyymmdd');
        
    case 'days'
        
        yearsPassed = floor(timepassed/365);
        extraDays = timepassed - yearsPassed*365;
        
        newDay = extraDays + start(3);
        if newDay > 365
            yearsPassed = yearsPassed + 1;
        end
        dayOfYear = mod(newDay,365);
        
        
        monthStart = [1 32 60 91 121 152 182 213 244 274 305 335];
        month = find(dayOfYear - monthStart >= 0, 1, 'last');
        day = dayOfYear - monthStart(month) + 1;
        day = floor(day);
        
        newDate = [start(1)+yearsPassed month day start(4:end)];
        str = datestr(newDate, 'yyyymmdd');
end

        
        
        
        