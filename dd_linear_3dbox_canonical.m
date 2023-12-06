function C = dd_linear_3dbox_canonical(c_f, U)

C = zeros(3,1);

[Gmat,Hmat] = get_hcoeffs(n,d);

h_coeffs = get_coeffs_3dbox(U);

temp_exps = homopoly(3,2);

I = 2*eye(3);
J = ones(3,3) - .5*I;
T = [2,3;1,3;1,2];
for i=1:3
    
    c_temp = 12*(4*Hmat(:,:,lex_index(I(i,:)))*c_f - abs(Hmat(:,:,lex_index(J(T(i,1),:)))*c_f) - abs(Hmat(:,:,lex_index(J(T(i,2),:)))*c_f));
    
    [val,idx] = max(-c_temp);
    
    C(i) = (max([val,0])/min(h_coeffs)) + 1;
    
end



