% Program: Socexp, WZJ, 20210721
% Input:  R   --- the payoff array, 
%        id   --- the number of the game 
% Output: the 5x5 payoff matrix
function retM5x5 = get_payoffmatrix_from_index(Ro,id)
	j=(id-1)*5; 
	retM5x5 = Ro(j+1:j+5,:);  
end
