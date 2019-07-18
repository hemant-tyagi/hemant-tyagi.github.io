function [dmin, match] = MatchLocations(test,t)

K = length(t);
dmin = zeros(K,1); % ith entry contains distance value for t_i 
match = zeros(K,1); %ith entry contains index of test which mathches with t_i

% Initialize distance matrix
distMat = zeros(K,K);

for i=1:K 
    for j = 1:K
    distMat(i,j) = min(abs(test(i) - t(j)), 1-abs(test(i) - t(j))); % wrap around distance
    end
end

for count = 1:K
   [curmin,row_ind] = min(distMat);
   [mindist, col_ind] = min(curmin); 
   row_ind = row_ind(col_ind);
   
   dmin(col_ind) = mindist;
   match(col_ind) = row_ind;   
   
   % this row and col will not take part in future calculations
   distMat(row_ind,:) = 10;
   distMat(:,col_ind) = 10;
end