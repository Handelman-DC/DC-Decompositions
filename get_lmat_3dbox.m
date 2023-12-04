function L = get_lmat_3dbox(u)

n1 = nchoosek(6,2);
n2 = nchoosek(3,2);

P = [0 0 0 1; u(3) 0 0 -1; 0 0 1 0; u(2) 0  -1 0; 0 1 0 0; u(1) -1 0 0]

tempL = zeros(n1,n2);

m_exps = lex_exps(3,2);
p_exps = lex_exps(6,2);

for i=1:n1;