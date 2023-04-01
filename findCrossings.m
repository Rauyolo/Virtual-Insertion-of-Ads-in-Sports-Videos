function [intersectionsForCalibration, intersectionsForMapping, caseNumber] = findCrossings(fileName, videoNum)

numImWithGoodIntersections = 1;

% create videoReader object
videoReader = VideoReader(fileName);
numFrames = videoReader.NumFrames;

% preallocation used for speed. truncated later
intersectionsForCalibration = cell(length(numFrames),1);
intersectionsForMapping = cell(length(numFrames),1);
caseNumber = zeros(1,length(numFrames));

f = waitbar(0,'1','Name',sprintf('Mapping points for video %d...',videoNum),...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

im_num = 0;
while hasFrame(videoReader)
    % update waitbar
    im_num = im_num + 1;
    waitbar(im_num/numFrames,f,sprintf('Processing image %d/%d',im_num,numFrames))
    
    % read current frame of the video
    original = readFrame(videoReader);

    % remove non-white objects
    outputImage = rgb2gray(original)>200;

    % remove noise (small objects)
    connectedComp = bwconncomp(outputImage,8);
    regProps = regionprops(connectedComp, 'Area');
    outputImage = ismember(labelmatrix(connectedComp),find([regProps.Area]>500))>0;
    
    % close olympic rings
    se = strel('disk',20);
    outputImage = imclose(outputImage,se);
    
    % find olympic rings by circularity and remove them
    connectedComp = bwconncomp(outputImage,8);
    regProps = regionprops(connectedComp,'Circularity');
    outputImage = ismember(labelmatrix(connectedComp),find([regProps.Circularity]<0.2))>0;
    
    % thin the result for ransac use later
    outputImage = bwmorph(outputImage,'thin', 'inf');   

    % obtain locations of all 1's
    [y, x] = find(outputImage);
    points = [x, y];
    
    % ransac expects 2*n matrix, transpose needed
    points = points.';

    % ransac is adapted from https://nl.mathworks.com/matlabcentral/fileexchange/30809-ransac-algorithm-with-example-of-finding-homography
    iterNum = 300;
    thDist = 2;
    thInlrRatio = .1;

    % slope and y-intercept values of lines. preallocation is used for
    % speed. truncated later
    k = zeros(1,10);
    b = zeros(1,10);
    angles = zeros(1,10);
    typeOfLine = zeros(1,10);

    % for plotting lines
    %figure(1);
    %imshow(outputImage); hold on;

    % find all slopes and y-intercepts of possible lines present in the
    % image
    count = 1;
    while 1
        % after first iteration remove previous line inliers from points
        if count ~= 1
           points(:,inliers) = []; 
        end

        % call ransac and to find k, b and points which are inliers
        [t,r,inliers] = ransac(points,iterNum,thDist,thInlrRatio);

        % no more lines detected => break
        if inliers == -1
            break;
        % smaller line was detected but it was found to be the net =>
        % remove it without adding it to the lines
        elseif t == -1
            continue;
        end
        
        % store values of k and b
        k(count) = -tan(t);
        b(count) = r/cos(t);
        angle = rad2deg(atan(k(count)));
        angles(count) = angle;

        % 1 is horizontal line, 0 is vertical
        if angle < 20 && angle >-20
            typeOfLine(count) = 1;
        else
            typeOfLine(count) = 0;
        end

        % plot the line
        %xVal = 1:1:size(original,2);
        %plot(points(1,:), points(2,:), 'b.'); hold on;
        %plot(xVal,k(count)*xVal+b(count),'y'); hold on;
        count = count + 1;
    end

    % removing unused preallocated values
    lastUsedIndex = find(k, 1, 'last');
    k = k(1:lastUsedIndex);
    b = b(1:lastUsedIndex);
    angles = angles(1:lastUsedIndex);
    typeOfLine = typeOfLine(1:lastUsedIndex);
        
    % find number of horiz and vertical lines
    numVertLines = nnz(~typeOfLine);
    numHorizLines = length(k) - numVertLines;

    % keep best 2 horizontal lines based on the number of points used to
    % construct them
    if numHorizLines > 2
        indexes = find(typeOfLine==1);
        k(indexes(3:end)) = [];
        b(indexes(3:end)) = [];
        angles(indexes(3:end)) = [];
        typeOfLine(indexes(3:end)) = [];
        numHorizLines = 2;
    end
    
    % finding intersections of lines
    allIntersections = zeros(1000,2);
    count = 1;
    for i=1:numHorizLines
        for j=numHorizLines+1:length(k)   
            % find x and y coordinates of the intersection
            xCoord = (b(i)-b(j))/(k(j)-k(i));
            yCoord = k(j)*xCoord + b(j);
            
            % if within image
            if xCoord <=length(outputImage) && yCoord <=width(outputImage) &&...
               xCoord >= 1 && yCoord >= 1
                % store intersections for this image
                allIntersections(count, :) = [xCoord, yCoord];
                count = count + 1;
                %plot(xCoord,yCoord,'x','LineWidth',20,'Color','green');
                %hold on;
            end
        end
    end

    % removing unused preallocated values
    allIntersections = allIntersections(any(allIntersections,2),:);

    %{
    % plotting intersections
    figure(1);
    imshow(outputImage);   hold on; 

    for i=1:length(allIntersections)
        plot(allIntersections(i,1),allIntersections(i,2),'x','LineWidth',20,'Color','green'); hold on;
    end

    title("Frame " + num2str(numImWithGoodIntersections));
    pause(0.1);
    %}

    % determine which world points are present in current image by case
    if length(allIntersections)==8
        caseNumber(numImWithGoodIntersections) = 1;
    elseif (length(allIntersections)==7 || length(allIntersections)==6) && angles(1) >= 0
        caseNumber(numImWithGoodIntersections) = 2;
    elseif (length(allIntersections)==7 || length(allIntersections)==6) && angles(1) < 0
        caseNumber(numImWithGoodIntersections) = 3;
    else
        caseNumber(numImWithGoodIntersections) = 4;
    end

    tempSort= sortrows(allIntersections,1);

    % if 7 points found, use instead 6 vertical pairs
    if caseNumber(numImWithGoodIntersections)==2 && length(allIntersections)==7
        tempSort(7,:) = [];
    elseif caseNumber(numImWithGoodIntersections)==3 && length(allIntersections)==7
        tempSort(1,:) = [];
    end

    % sorting points so that they map to the world points
    intersectionSorted = [];
    for i=1:2:length(allIntersections)-1
        temp = sortrows(tempSort(i:i+1,:),2);
        intersectionSorted = cat(1,intersectionSorted,temp);
    end
    
    % add values of NaN to world points not in calibration images (see
    % estimateCameraParameters)
    useForCalibration = 1;
    if caseNumber(numImWithGoodIntersections)==1
        intersectionSortedCalibration = intersectionSorted;
    elseif caseNumber(numImWithGoodIntersections)==2
        intersectionSortedCalibration = [[NaN, NaN]; [NaN, NaN]; intersectionSorted];
    elseif caseNumber(numImWithGoodIntersections)==3
        intersectionSortedCalibration = [intersectionSorted; [NaN, NaN]; [NaN, NaN]];
    elseif caseNumber(numImWithGoodIntersections)==4
        useForCalibration = 0;
    end

    % if not 6 or 8 points, don't use this point for calibration
    if useForCalibration
        [intersectionsForCalibration{numImWithGoodIntersections}(:,:)] = intersectionSortedCalibration;
    end

    % store found intersections for mapping
    [intersectionsForMapping{numImWithGoodIntersections}(:,:)] = intersectionSorted;

    % count these intersections as good for indexing
    numImWithGoodIntersections = numImWithGoodIntersections + 1;
end

% removing empty cells
intersectionsForCalibration = intersectionsForCalibration(~cellfun('isempty',intersectionsForCalibration));
intersectionsForMapping = intersectionsForMapping(~cellfun('isempty',intersectionsForMapping));
caseNumber = caseNumber(1:find(caseNumber, 1, 'last'));

% delete waitbar
delete(f)

% close open figures
close all;
end