function steerFiltHyperVolume = makeSteerBasis3D(filt, steerDirectionArray)

    nfilt = size(steerDirectionArray, 1);
    [nrows, ncols, nframes] = size(filt);
    steerFiltHyperVolume = zeros(nrows, ncols, nframes, nfilt);
    for ifilt = 1:nfilt
        steerDir = steerDirectionArray(ifilt,:);
        steerFilt = rotateFilt3D(filt, steerDir);
        steerFiltHyperVolume(:,:,:,ifilt) = steerFilt;
    end
end

