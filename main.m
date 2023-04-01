close all
clear variables

% read logo image and add a frame around it
logoImage = imread("rings.png");
logoImage = padarray(logoImage,[20 20],0,'both');
logoImage = padarray(logoImage,[5 5],255,'both');
PSF = fspecial('gaussian',120,20);
logoImage = edgetaper(logoImage,PSF);
[lengthLogo, widthLogo, ~] = size(logoImage);
logoPoints = [1, 1; widthLogo, 1; widthLogo, lengthLogo; 1, lengthLogo];

% manually found intersection points in video 1
load("actualCrossings.mat")

for videoNum = 1:4
    fileName = strcat('video', num2str(videoNum), '.mp4');
    videoReader = VideoReader(fileName);
    numFrames = videoReader.NumFrames;

    % find intersections for calibration and mapping
    [intersectionsForCalibration, intersectionsForMapping, caseNumber] = findCrossings(fileName, videoNum);

    figure;
    vanishingPointError = findVanishingPointError(intersectionsForMapping);
    bar(vanishingPointError);
    xlabel('Frame number')
    ylabel('Mean pixel distance to vanishing point')
    title("Vanishing point error for video " + num2str(videoNum))

    % locations of intersections (in meters)
    worldPoints = [0.0,0.0; 0.0,9.0; 6.0,0.0; 6.0,9.0; 12.0,0.0; 12.0,9.0; 18.0,0.0; 18.0,9.0];
    
    % do camera calibration only with video 1
    if videoNum == 1
        % calibrate camera
        [cameraParams, validImageIndices, estimationErrors] = ...
            estimateCameraParameters(cat(3,intersectionsForCalibration{:}), worldPoints, ...
            'ImageSize', [videoReader.height, videoReader.width], 'WorldUnits', 'm');

        % intrinsic camera parameters
        intrinsics = cameraParams.Intrinsics;
    end
    
    mappedImPts = zeros(4,2,numFrames);
    
    % logo locations in world coordinates
    placementWorldPoints(:,:,1) = [0.5,-4.0,0.0; 5.5,-4.0,0.0; 5.5,-1.0,0.0; 0.5,-1.0,0.0];
    %placementWorldPoints(:,:,2) = [12.5, -1.245,-0.3; 15.5,-1.245,-0.3; 15.5,-1.0,0.0; 12.5,-1.0,0.0];
    placementWorldPoints(:,:,2) = [19.5,2.0,-0.3; 19.5,7.0,-0.3; 18.5,7.0,0.0; 18.5,2.0,0.0];
    
    % for extrinsics
    for i=1:length(intersectionsForMapping)
        % find appropriate worldPoints based on the case
        % if less than 6 points are detected (case 4), hard to know which  
        % points were found, so instead the previous extrinsic values are  
        % used again, assuming that the camera did not move much in new frame.
        notEnoughIntersections = 0;
        switch caseNumber(i)
        case 1
            worldPoints = [0.0,0.0; 0.0,9.0; 6.0,0.0; 6.0,9.0; 12.0,0.0; 12.0,9.0; 18.0,0.0; 18.0,9.0];
        case 2
            worldPoints = [0.0,0.0; 0.0,9.0; 6.0,0.0; 6.0,9.0; 12.0,0.0; 12.0,9.0];
        case 3
            worldPoints = [6.0,0.0; 6.0,9.0; 12.0,0.0; 12.0,9.0; 18.0,0.0; 18.0,9.0];
        case 4
            notEnoughIntersections = 1;
        end
    
        % if normal case, estimate extrinsics
        if ~notEnoughIntersections && vanishingPointError(i)~=-1
            camExtrinsics = estimateExtrinsics(cell2mat(intersectionsForMapping(i)),worldPoints,intrinsics);
        end
    
        for logoNum = 1:size(placementWorldPoints,3)
            % mapping destination world coordinates to image coordinates
            mappedImPts(:,:,i,logoNum) = world2img(placementWorldPoints(:,:,logoNum),camExtrinsics,intrinsics);
        end
    end
    
    % get corresponding output file for this video
    outputFile = strcat('output',num2str(videoNum));

    % create the video object
    video = VideoWriter(outputFile, 'MPEG-4'); 

    % open the file for writing  
    open(video); 
    
    % create waitbar
    f = waitbar(0,'1','Name',sprintf('Placing logo for video %d...',videoNum),...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

    im_num = 0;
    while hasFrame(videoReader)
        % update waitbar
        im_num = im_num + 1;
        waitbar(im_num/numFrames,f,sprintf('Processing image %d/%d',im_num,numFrames))

        % read current video frame
        courtImage = readFrame(videoReader);
        
        for logoNum = 1:size(mappedImPts,4)
            % load mapped location of logo for this image
            mappedPoints = mappedImPts(:, :, im_num, logoNum);
            
            % perform homography and find H
            tform = estgeotform2d(mappedPoints, logoPoints, 'projective');
            H = tform.A;
            
            % find coordinates of the pixels inside the region of mappedImPts
            pixelsInsideLogo = poly2mask(mappedPoints(:, 1), mappedPoints(:, 2), size(courtImage,1), size(courtImage,2));
            [xCoord, yCoord] = find(pixelsInsideLogo==1);
            pixelsInsideLogo = [yCoord, xCoord, ones(length(xCoord),1)];
            
            % obtain warped points locations of the world image
            product=H*pixelsInsideLogo';
            warpedLocations = [product(1,:)./product(3,:); product(2,:)./product(3,:)]';
            
            % round up the values of warpedLocations
            warpedLocations = ceil(warpedLocations);
            
            % fix rounding errors caused by Matlab
            if any(min(warpedLocations(:)) < 1) ||...
                    (any(max(warpedLocations(:, 1)) > size(logoImage,2)) ||...
                    any(max(warpedLocations(:, 2)) > size(logoImage,1)))
                warpedLocations(warpedLocations(:) < 1) = 1;
                warpedLocations(warpedLocations(:,1) > size(logoImage,2), 1) = size(logoImage,2);
                warpedLocations(warpedLocations(:,2) > size(logoImage,1), 2) = size(logoImage,1);
            end

            % if less matlab warping function has some bugs
            if length(warpedLocations)<50
                continue;
            end

            % find world limits for the logo
            qq = [size(logoImage,1), 0;
                  size(logoImage,1), size(logoImage,2);
                  0,                 size(logoImage,2);
                  0,                 0];

            qq_xWorldLimits = [min(qq(:,1)), max(qq(:,1))];
            qq_yWorldLimits = [min(qq(:,2)), max(qq(:,2))];
            ad_img_ref2 = imref2d(size(logoImage), qq_xWorldLimits*2, qq_yWorldLimits/2);

            % find tform for warping
            tform1 = estgeotform2d(warpedLocations, pixelsInsideLogo(:,1:2), 'projective');

            % warp the logo to destination location
            warpedLogo = imwarp(logoImage, ad_img_ref2, tform1, 'OutputView', imref2d(size(courtImage)));

            % apply thresholding on the HSV to get the mask
            [courtMask, ~] = createMask(courtImage);

            % apply mask to image
            masked_warped_logo = im2uint8(im2double(warpedLogo) .* ~courtMask);
            idx = (sum(masked_warped_logo, 3) ~= 0);
            idx = logical(idx.*ones([size(idx),3]));
            projected_img = courtImage;
            projected_img(idx) = masked_warped_logo(idx);
            courtImage = projected_img;
        end

        % add red circles at locations of detected intersections
        switch caseNumber(im_num)
        case 1
            label = ["(0.0,0.0)"; "(0.0,9.0)"; "(6.0,0.0)"; "(6.0,9.0)"; "(12.0,0.0)"; "(12.0,9.0)"; "(18.0,0.0)"; "(18.0,9.0)"];
        case 2
            label = ["(0.0,0.0)"; "(0.0,9.0)"; "(6.0,0.0)"; "(6.0,9.0)"; "(12.0,0.0)"; "(12.0,9.0)"];
        case 3
            label = ["(6.0,0.0)"; "(6.0,9.0)"; "(12.0,0.0)"; "(12.0,9.0)"; "(18.0,0.0)"; "(18.0,9.0)"];
        case 4
            label = "Error";
        end
    
        for i=1:length(cell2mat(intersectionsForMapping(im_num)))
            temp = cell2mat(intersectionsForMapping(im_num));
            position = [temp(i,1),temp(i,2),10];

            if length(label)==1
                currentLabel = label;
            else
                currentLabel = label(i);
            end
            
            courtImage = insertObjectAnnotation(courtImage,"circle",position,currentLabel,LineWidth=3,Color="red",TextColor="black");
        end
        
        % write the image to file
        writeVideo(video,courtImage);
        imshow(courtImage);

        % calculate reprojection error for video 1 only
        if videoNum == 1
            cumulative = 0;
            offset = 0;
    
            switch caseNumber(im_num)
            case 1
                numPoints = 4;
            case 2
                numPoints = 3;
            case 3
                numPoints = 3;
            case 4
                numPoints = 1;
            end
    
            for point = 1:4
                if actualCrossings{point, im_num} == 0
                    offset = offset + 1;
                    continue
                else
                    pointsAt = actualCrossings{point, im_num};
                    pointsAt2 = actualCrossings{point + 4, im_num};
    
                    difference = pointsAt - intersectionsForMapping{1, im_num}(2*(point-offset-1) + 1, :);
                    difference2 = pointsAt2 - intersectionsForMapping{1, im_num}(2*(point-offset-1) + 2, :);
    
                    cumulative = cumulative + sqrt(difference(1,1)^2 + difference(1, 2)^2) + sqrt(difference2(1,1)^2 + difference2(1, 2)^2);
                end
            end
            cumulative = cumulative / (2 * numPoints);
            diffVector(im_num, 1) = cumulative;
        end
        %
    end
    
    % close the file
    close(video); 
    %close;
    delete(f);

    % plot reprojection error for video 1
    if videoNum == 1
        figure;
        plot(1:1:numFrames, diffVector);
        xlabel('Video frame number');
        ylabel('Average point error per frame [px]')
        title('Reprojection error (detected line crossing errors)')
    end
end