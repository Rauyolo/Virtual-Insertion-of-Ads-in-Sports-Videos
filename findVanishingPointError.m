function vanishingPointError = findVanishingPointError(intersectionsForMapping)
    vanishPtsEst = cell(length(intersectionsForMapping),1);
    vanishingPointError = zeros(length(intersectionsForMapping),1);

    % loop through intersections of every frame
    for i=1:length(intersectionsForMapping)
        intersections = cell2mat(intersectionsForMapping(i));
        k = zeros(1, length(intersections)/2);
        b = zeros(1, length(intersections)/2);

        % find equations of vertical lines 
        for j=1:2:length(intersections)
            index = nnz(k)+1;
            k(index) = (intersections(j+1,2)-intersections(j,2))/(intersections(j+1,1)-intersections(j,1));
            b(index) = intersections(j,2)-k(index)*intersections(j,1);
        end

        vanishPtsEst{i}(:,:) = zeros(length(k),2);
        count = 1;

        % find vanishing points (intersections of vertical lines)
        for ii=1:length(k)
            for j=ii+1:length(k)
                xCoord = (b(ii)-b(j))/(k(j)-k(ii));
                yCoord = k(j)*xCoord + b(j);
                vanishPtsEst{i}(count,:) = [xCoord, yCoord];
                count = count + 1;
            end
        end

        % if less than 3 lines => error frame (caseNumber=4)
        if length(vanishPtsEst{i}) < 3
            vanishingPointError(i) = -1;
            continue;
        else
            % find mean of intersections of lines
            meanOfPoints = [mean(vanishPtsEst{i}(:,1)), mean(vanishPtsEst{i}(:,2))];
        end
        
        % find smallest distance between all vertical lines and mean
        xCoords = (meanOfPoints(1)+k.*(meanOfPoints(2)-b))./(k.^2+1);
        yCoords = k.*xCoords+b;
        disLinesToMean = sqrt((xCoords-meanOfPoints(1)).^2+(yCoords-meanOfPoints(2)).^2);
        error = mean(disLinesToMean);

        % if error is large, probably points got shuffled around => error
        if error > 100
           vanishingPointError(i) = -1;
        else
           vanishingPointError(i) = error;
        end
    end
end