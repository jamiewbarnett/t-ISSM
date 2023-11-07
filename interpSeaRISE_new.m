function [dataout] = interpSeaRISE(X,Y,string,varargin),
%INTERPSEARISE - interpolate SeaRISE data
%
%   Available data:
%      1.  sealeveltimes
%      2.  dhdt
%      3.  surfvelmag
%      4.  balvelmag
%      5.  oisotopestimes
%      6.  bheatflx
%      7.  presprcp
%      8.  sealevel_time_series
%      9.  usrf
%      10. mapping
%      11. surfvely
%      12. surfvelx
%      13. topg
%      14. landcover
%      15. temp_time_series
%      16. thk
%      17. time
%      18. oisotopes_time_series
%      19. runoff
%      20. smb
%      21. airtemp2m
%      22. surftemp
%
%   Usage:
%      [dataout] = interpBedmap2(X,Y,string)

verbose=0;

if nargin==3,
	hemisphere = +1;
else
	hemisphere = varargin{1};
end


%read data

if hemisphere==1,
	searisenc='./Model_Data/Greenland_5km_dev1.2.nc';
elseif hemisphere==-1,
	searisenc='/home/ModelData/SeaRISE/Antarctica_5km_dev1.0.nc';
end
	

%convert coordinates to SeaRISE projection
if verbose, disp('   -- SeaRISE: converting coordinates'); end
if hemisphere==1,
	[LAT,  LON  ] = xy2ll(double(X(:)),double(Y(:)),+1,45,70);
	[xproj,yproj] = ll2xy(LAT,LON  ,+1,39,71);
elseif hemisphere==-1,
	xproj=X; yproj=Y;
end

if verbose, disp('   -- SeaRISE: loading coordinates'); end
xdata = double(ncread(searisenc,'x1'));%*1000;
ydata = double(ncread(searisenc,'y1'));%*1000;

if verbose, disp(['   -- SeaRISE: loading ' string]); end
data  = double(ncread(searisenc,string))';

if verbose, disp(['   -- SeaRISE: interpolating ' string]); end
if strcmpi(string,'LandMask');
	dataout = InterpFromGrid(xdata,ydata,data,xproj,yproj,'nearest');
else
	dataout = InterpFromGrid(xdata,ydata,data,xproj,yproj);
end
dataout = reshape(dataout,size(X,1),size(X,2));
