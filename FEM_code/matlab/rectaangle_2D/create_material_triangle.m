function c = create_material_triangle( ...
    E,miu,thickness)
    c = thickness*E/(1-miu^2)*[1,miu,0;miu,1,0;0,0,(1-miu)/2];
end