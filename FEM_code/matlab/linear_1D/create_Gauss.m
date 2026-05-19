function [Gauss_weights,Gauss_points] = create_Gauss(type_dimesion,number_order)
    n = round((number_order+1)/2);
    if n==1
        Gauss_weights_1D = 2;
        Gauss_points_1D = 0;
    elseif n==2
        Gauss_weights_1D = [1,1];
        Gauss_points_1D = [-sqrt(1/3),sqrt(1/3)]; 
    end

    if type_dimesion==1
        Gauss_weights = Gauss_weights_1D;
        Gauss_points = Gauss_points_1D;
    end
end