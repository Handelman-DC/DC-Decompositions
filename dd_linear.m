function C = dd_linear(c_f, U)

C = zeros(3,1);

[Gmat,Hmat] = get_hcoeffs(n,d);

h_coeffs = get_coeffs_poly3dbox(U);

[L,Linv] = get_Lmat_3dbox(U);

temp_exps = homopoly(6,2);

I = 2*eye(3);
J = ones(3,3) - .5*I;
T = [2,3;1,3;1,2];
for i=1:6
    
    c_temp = 12*(Linv*Hmat(:,:,lex_index(I(i,:)))*c_f - abs(Linv*Hmat(:,:,lex_index(J(T(i,1),:)))*c_f) - abs(Linv*Hmat(:,:,lex_index(J(T(i,2),:)))*c_f));
    
    [val,idx] = max(-c_temp);
    
    C(i) = (max([val,0])/min(h_coeffs)) + 1;
    
end



