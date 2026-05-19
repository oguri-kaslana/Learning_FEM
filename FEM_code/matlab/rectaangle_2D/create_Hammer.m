function [Hammer_points,Hammer_weights] = create_Hammer(number_dimesion,number_order)
    if number_dimesion ~= 2
        error('当前 create_Hammer 只用于二维三角形单元');
    end

    if number_order == 1

        Hammer_points = [
            1/3, 1/3
        ];

        Hammer_weights = [
            1/2
        ];

    elseif number_order == 2

        Hammer_points = [
            1/6, 1/6;
            2/3, 1/6;
            1/6, 2/3
        ];

        Hammer_weights = [
            1/6;
            1/6;
            1/6
        ];

    elseif number_order == 3

        Hammer_points = [
            1/3, 1/3;
            1/5, 1/5;
            3/5, 1/5;
            1/5, 3/5
        ];

        Hammer_weights = [
            -27/96;
             25/96;
             25/96;
             25/96
        ];

    else

        error('暂时只支持 number_order = 1, 2, 3');

    end

end