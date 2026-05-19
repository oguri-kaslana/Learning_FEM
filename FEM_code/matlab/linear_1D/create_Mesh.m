function [point,topolocy,number_elements]=create_Mesh( ...
    number_dimesion,...
    left,right,step)
    if number_dimesion ==1
        point = (left:step:right).';
        number_elements = round((right-left)/step);
        topolocy = zeros(number_elements,2);
        for i=1:number_elements
            topolocy(i,:) = [i,i+1];
        end
    end
end