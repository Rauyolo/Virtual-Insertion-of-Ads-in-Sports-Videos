function index = randIndex(maxIndex,len)
%INDEX = RANDINDEX(MAXINDEX,LEN)
%   randomly, non-repeatedly select LEN integers from 1:MAXINDEX

%{ 
Adapted from:
Ke Yan (2022). RANSAC algorithm with example of finding homog-
raphy (https://www.mathworks.com/matlabcentral/fileexchange/30809-
ransac-algorithm-with-example-of-finding-homography), MATLAB
Central File Exchange. Retrieved November 4, 2022
%}

if len > maxIndex
	index = [];
	return
end

index = zeros(1,len);
available = 1:maxIndex;
rs = ceil(rand(1,len).*(maxIndex:-1:maxIndex-len+1));
for p = 1:len
	while rs(p) == 0
		rs(p) = ceil(rand(1)*(maxIndex-p+1));
	end
	index(p) = available(rs(p));
	available(rs(p)) = [];
end