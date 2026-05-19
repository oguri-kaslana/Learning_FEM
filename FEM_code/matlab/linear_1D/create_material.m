function c = create_material( ...
    number_dimesion, ...
    E,A,miu)
    if number_dimesion==1
        c = E*A;
    end
end