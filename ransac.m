function [ theta,rho,max_inliers ] = ransac( pts,iterNum,thDist,thInlrRatio )
%RANSAC Use RANdom SAmple Consensus to fit a line
%	RESCOEF = RANSAC(PTS,ITERNUM,THDIST,THINLRRATIO) PTS is 2*n matrix including 
%	n points, ITERNUM is the number of iteration, THDIST is the inlier 
%	distance threshold and ROUND(THINLRRATIO*SIZE(PTS,2)) is the inlier number threshold. The final 
%	fitted line is RHO = sin(THETA)*x+cos(THETA)*y.
%	Yan Ke @ THUEE, xjed09@gmail.com

%{ 
Adapted from:
Ke Yan (2022). RANSAC algorithm with example of finding homog-
raphy (https://www.mathworks.com/matlabcentral/fileexchange/30809-
ransac-algorithm-with-example-of-finding-homography), MATLAB
Central File Exchange. Retrieved November 4, 2022

Edited by Raul Ismayilov
%}

sampleNum = 2;
ptNum = size(pts,2);
thInlr = round(thInlrRatio*ptNum);
inlrNum = zeros(1,iterNum);
theta1 = zeros(1,iterNum);
rho1 = zeros(1,iterNum);
max_inliers = [];
for p = 1:iterNum
	% 1. fit using 2 random points
	sampleIdx = randIndex(ptNum,sampleNum);
	ptSample = pts(:,sampleIdx);
	d = ptSample(:,2)-ptSample(:,1);
	d = d/norm(d); % direction vector of the line
	
	% 2. count the inliers, if more than thInlr, refit; else iterate
	n = [-d(2),d(1)]; % unit normal vector of the line
	dist1 = n*(pts-repmat(ptSample(:,1),1,ptNum));
	inlier1 = find(abs(dist1) < thDist);
    if length(max_inliers) < length(inlier1)
        max_inliers = inlier1;
    end
	inlrNum(p) = length(inlier1);
	if length(inlier1) < thInlr, continue; end
	ev = pca(pts(:,inlier1)');
	d1 = ev(:,1);
	theta1(p) = -atan2(d1(2),d1(1)); % save the coefs
	rho1(p) = [-d1(2),d1(1)]*mean(pts(:,inlier1),2);
end

% 3. choose the coef with the most inliers
if length(max_inliers) > 150
    [~,idx] = max(inlrNum);
    theta_temp = theta1(idx);
    angle = rad2deg(atan(-tan(theta_temp)));
    if (angle < 0 && angle > -85) || (angle > 0 && angle < 85)
        theta = theta_temp;
        rho = rho1(idx);
    else
        theta = -1;
        rho = -1;
    end
elseif length(max_inliers) > 80
    [~,idx] = max(inlrNum);
    theta_temp = theta1(idx);
    angle = rad2deg(atan(-tan(theta_temp)));
    if (angle < -20 && angle > -80) || (angle > 20 && angle < 80)
        theta = theta_temp;
        rho = rho1(idx);
    else
        theta = -1;
        rho = -1;
    end
else
    theta = -1;
    rho = -1;
    max_inliers = -1;
end

end