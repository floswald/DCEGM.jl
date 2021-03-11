% Handy function to display time in human readable form
% Written by Fedor Iskhakov, University of New South Wales, 2015

function str = ht(x)
% This function prints elapsed time in human readable format
% x= time in seconds 
% Created by: Fedor Iskhakov

names={'s','m','h',' day',' week',' year'};
plural={'','','','s','s','s'};
dividers=[1,60,60*60,24*60*60,7*24*60*60,52*7*24*60*60];

str='';

for i=numel(names):-1:1
	if dividers(i)>1
		dig=floor(x/dividers(i));
		digstr=sprintf('%1.0f',dig);		
	else
		dig=x;
		digstr=sprintf('%1.3f',dig);		
	end
	if dig>1
		str=sprintf('%s %s%s%s',str,digstr,names{i},plural{i});
	elseif dig<=1 && dig>0
		str=sprintf('%s %s%s',str,digstr,names{i});
	end
	x=x-dig*dividers(i);
end