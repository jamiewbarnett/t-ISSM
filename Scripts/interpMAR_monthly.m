function output = interpMAR_monthly(X,Y,string, time, filename)
%%%%%% %Interpolate from MAR (Fettweis et al 2017)

%Available outputs
	% -LAT
	% -LON
	% -ME					(Melt)
	% -MSK_BAM01		(Mask Bamber et al 2001 5x5 km)
	% -MSK_BAM13		(Mask Bamber et al 2013 1x1 km)
	% -MSK_MAR			(MAR 10x10 km Ice Mask)
	% -RF					(Rainfall)
	% -RU					(Runoff)
	% -SF					(Snowfall)
	% -SMB				(Surface Mass Balance (without corrections)
	% -SMBCORR			(Surface Mass Balance (with corrections)
	% -SRF_BAM01		(Bamber et al 2001 5x5km Surface height)
	% -SRF_BAM13		(Bamber et al 2013 1x1 km Surface height)
	% -SRF_MAR			(MAR 10x10 km Surface height)
	% -ST					(Surface temperature)
	% -ST					(Surface temperature (with corrections)
	% -SU					(Sublimation/Evaporation)
	% -TIME				(Time)
	% -X					(x)
	% -Y					(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isverbose=0;

%Define data source
%ncpath='/home/fholmes/Documents/issm_scripts/ModelData/MARv3.5.2-10km-yearly-ERA-Interim-1979-2014.nc';
ncpath= filename;

% Build grid from Bamber2001
xdata = -760000:20000:680000;
ydata = -3330000:20000:-640000;
%ydata = -1180000:20000:1500000; 

%Convert to MAR projections
if isverbose, disp('   -- MAR: converting coordinates'); end
[LAT,LON] = xy2ll(double(X(:)),double(Y(:)),+1,45,70); %convert model mesh xy to Bamber's projection
[xMAR,yMAR] = ll2xy(LAT,LON,+1,39,71); %convert from lat/long Bamber's to MAR projection xy

if isverbose, disp('   -- MAR: loading data'); end
data  = double(ncread(ncpath,string));

%Define what year to load data for 1979-2014
starttime=1;
year=time-starttime+1;

%Get data from given year
data2=squeeze(data(:,:,year));

if isverbose, disp('   -- MAR: interpolating data to grid'); end
output = InterpFromGrid(xdata,ydata,data2',xMAR,yMAR);
output = reshape(output,size(X,1),size(X,2));
