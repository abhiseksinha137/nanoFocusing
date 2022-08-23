function []=displayLoop(current, total)
disp([num2str(current/total*100, '%2.2f'), '% Complete...']);