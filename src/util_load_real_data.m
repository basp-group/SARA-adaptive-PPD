function [visibilities, ucoorRadians, vcoorRadians,wcoorUnitWavelength, nWeights,time,pixelSize] = util_load_real_data(visibilityFileName, param)
%% cst
speedOfLight = 299792458;
%% read param
obsFreq = param.obsFreq;
pixelSize = param.pixelSize; %in arcsec
imageResolutionStr = param.imageResolutionStr;
%% load data
try load(visibilityFileName,'y_I','flag','weights','uvw','time') ;% flag ; sigmas ; uvw; y_I; weights; time
catch,   load(visibilityFileName,'y_I','flag','weights','uvw');
end
%
flag = double(flag(:));
y_I  = double(y_I(:));
weights = double(weights(:)); %1/(sigma^2)
%
ind2keep =(weights>0).*(flag==0).*(abs(y_I)>0);
%
if param.isWeightsFixed2Sigma
    weights = 1./(weights(ind2keep>0));
else
    weights = sqrt(weights(ind2keep>0));
end
%
y_I  = double(y_I(ind2keep>0));


% uvw
ucoor =  uvw(:,1);
ucoor = ucoor(ind2keep>0);

vcoor = -uvw(:,2);
vcoor = vcoor(ind2keep>0);

w = uvw(:,3);
wcoor = w(ind2keep>0);

if ~exist('time')
    time = ones(size(ucoor));
    fprintf("\nWarning: Time information needed for data blocking!\n")
end
time = double(time(ind2keep>0));
%% whitening the data: natural weighting
visibilities = y_I(:).* weights(:);
nWeights = weights;

%% Setting image resolution
maxProjBaseline = max(sqrt(ucoor.^2+vcoor.^2)); % maximum baseline
obsWavelength = speedOfLight/obsFreq;
if ~isempty(pixelSize)
    imagedFourierBandwidthScale =(obsWavelength*(180*60*60)/pi)/(2*(maxProjBaseline)*pixelSize); % resolution factor
else
    switch imageResolutionStr
        case 'nominal'
            imagedFourierBandwidthScale = 1.5;
            pixelSize = (obsWavelength*(180*60*60)/pi)/(imagedFourierBandwidthScale*2*(maxProjBaseline));
            disp(['INFO: pixelsize: ', num2str(pixelSize), ' asec ']);
        otherwise
            error('User should provide pixelsize')
    end
end
disp(['INFO: Imaging up to ',num2str(imagedFourierBandwidthScale), 'x the nominal resolution'] );
vcoorRadians = vcoor *pi/(maxProjBaseline * imagedFourierBandwidthScale); % scaling of u,v points for NUFFT
ucoorRadians = ucoor *pi/(maxProjBaseline * imagedFourierBandwidthScale); % scaling of u,v points for NUFFT
wcoorUnitWavelength = wcoor/obsWavelength;
