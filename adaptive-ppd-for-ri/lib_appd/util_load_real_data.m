function [y, uw, vw, sigma_noise, nW,time] = util_load_real_data(visibility_file_name, param)

load(visibility_file_name) % flag ; sigmas ; uvw; y_I; weights; time
flag = flag(:);
weights = double(weights(:));
ind=(weights>0).*(flag==0);
weights = sqrt(weights(ind>0));
y_I=double(y_I(ind>0));

% sigma_noise = 1.253*mad(y_V,0); %MAD estimator
sigma_noise = sqrt(2);

u =  uvw(:,1);
u =  u(ind>0);
v = -uvw(:,2);
v = v(ind>0);

%ant1 =double(ant1(ind>0));
%ant2 =double(ant2(ind>0));
%ant1=ant1(:);
%ant2=ant2(:);
if ~exist('time')
    time = ones(size(u));
else
    time =double(time(ind>0));
end

% whitening the data: natural weighting
y = y_I(:) .* weights(:);

% Setting imaging params
bmaxProj = max(sqrt(u.^2+v.^2)); % maximum baseline

%% resolution 
freq = param.freq;
cellsize = param.pixel_size; %in arcsec
speed_light = 299792458;
wavelength = speed_light/freq;
dl =(1.*wavelength*(180*60*60)/pi)/(2*(bmaxProj)*cellsize); % resolution factor
disp(['INFO: Imaging up to ',num2str(dl), 'x the nominal resolution'] );
disp ''
vw = v *pi/(bmaxProj * dl); % scaling of u,v points for NUFFT
uw = u *pi/(bmaxProj * dl); % scaling of u,v points for NUFFT
nW = weights;

