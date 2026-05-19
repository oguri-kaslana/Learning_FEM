function results = cauclate_results( ...
    number_dimesion,points,topolocy,...
    k_mod,f,c,k_global)
    %cauclate u
    u = cauclate_U(k_mod,f);

    %cauclate strain/stress
    [strain,stress] = cauclate_strain_stress(number_dimesion,points,topolocy,u,c);

    %cauclate RF
    RF = cauclate_RF(k_global,u,f);

    results = [u,RF,strain,stress];

end