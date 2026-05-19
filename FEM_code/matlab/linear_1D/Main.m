%create mesh
clc;
number_dimesion = 1;
[points,topolocy,number_elements] = create_Mesh(number_dimesion,0,2,1);

%create material
c = create_material(number_dimesion,100,1);

%translate
Jacobi = create_Jacobi(number_dimesion,points,topolocy);

%create B
B = create_B(number_dimesion);

%create K_local
k_local = create_K_local(number_dimesion,Jacobi,c);

%create K_global
k_global = create_K_global(number_dimesion,points,topolocy,Jacobi,c);

%create_boundary
q_fun = @(x) 10+2*x;
qb = Boundary_Q(1,2,points,topolocy,q_fun,1);
[k_mod,f] = Boundary_Disp(number_dimesion,1,0,k_global,qb);

%cauclate results(u,RF,strain,stress)
results = cauclate_results(number_dimesion,points,topolocy,k_mod,f,c,k_global);
disp(results);



