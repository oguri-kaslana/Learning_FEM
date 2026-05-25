function [Gp,Gw] = calculate_Gauss(number_Gp)
    if number_Gp == 1
        Gp = 0;
        Gw = 3;
    elseif number_Gp == 2
        Gp = [-sqrt(1/3);sqrt(1/3)];
        Gw = [1;1];
    elseif number_Gp == 3
        Gp = [-sqrt(3/5);0;sqrt(3/5)];
        Gw = [5/9;8/9;5/9];
    end

end