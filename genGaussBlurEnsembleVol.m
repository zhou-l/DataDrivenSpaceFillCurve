function genGaussBlurEnsembleVol(filename)
close all;
    headerInfo = nhdr_nrrd_read(filename, true);
    [filepath, name, ext] = fileparts(filename);
%     ltfilename = strcat(filepath, '/LT', name, '.csv');
%     vofilename = strcat(filepath, '/VO', name, '.csv');
    headerInfo = nhdr_nrrd_read(filename, true);
%     vol = single(headerInfo.data);
      vol = headerInfo.data;
    dim = size(vol);
    figure
    montage(vol);

    siz = dim;
    volx = squeeze(vol);
    for s = 0.5:0.5:5
        sigma = s;
        i = int16(s / 0.5);
        volSmooth = imgaussfilt3(volx, sigma);

%         figure
        montage(reshape(volSmooth,siz(1),siz(2),1,siz(3)))
        title('Gaussian filtered image volume');
        filtName = sprintf('%sSig%d.nhdr', name, i);
        nrrdWriter(filtName, volSmooth, [1 1 1], [0 0 0], 'raw');
    end

end