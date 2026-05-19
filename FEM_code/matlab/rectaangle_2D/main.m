clc;
clear;

[points,topolocy,nodes] = create_mesh_triangle(0,3,0,2);

[jacobi,det_jacobi] = create_jacobi_triangle(topolocy,points);

B = create_B_triangle(points,topolocy,jacobi);

c = create_material_triangle(1000,0.25,1);

k_local = create_K_local(topolocy,points,B,c,det_jacobi);

k_global = create_K_global(topolocy,points,k_local);

%f = Boundary_f(points,topolocy,[4,8,12],[100,0;100,0;100,0]);

%[k_mod,f_mod] = Boundary_u(points,topolocy,k_global,f,[1,5,9],[0,0;0,0;0,0]);

q = Boundary_q(points,topolocy,1,[12,8,4],[0,100;0,100]);

[k_mod,f_mod] = Boundary_u(points,topolocy,k_global,q,[1,5,9],[0,0;0,0;0,0]);

u = cauclate_U(k_mod,f_mod);

disp('节点位移 u = ');
disp(u);