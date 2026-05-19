function [jacobi,det_jacobi] = create_jacobi_triangle(topolocy,points)
    x1 = points(topolocy(:,1),1);
    x2 = points(topolocy(:,2),1);
    x3 = points(topolocy(:,3),1);
    y1 = points(topolocy(:,1),2);
    y2 = points(topolocy(:,2),2);
    y3 = points(topolocy(:,3),2);
    number_elements = size(topolocy,1);
    jacobi = zeros(2,2,number_elements);
    det_jacobi = zeros(number_elements,1);
    for i = 1: number_elements
        jacobi(:,:,i) = [
            x1(i)-x3(i),y1(i)-y3(i);
            x2(i)-x3(i),y2(i)-y3(i)];
        det_jacobi(i) = det(jacobi(:,:,i));
    end

end