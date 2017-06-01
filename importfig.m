function [x,y] = importfig
D=get(gca,'Children'); %get the handle of the line object
X=get(D,'XData'); %get the x data
Y=get(D,'YData'); %get the y data
x = [X{3} X{2} X{1}];
y = [Y{3} Y{2} Y{1}];
end

