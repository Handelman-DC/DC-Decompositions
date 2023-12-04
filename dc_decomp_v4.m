
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Local DC Decomposition of nonconvex polynomials using Handelman's Theorem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The main objective of this script is to accept a polynomial P(x), and
%returna two convex polynomials gamma(x) and nu(x) such that P = gamma-nu
%and both gamma and nu are convex. We do so using LP/SDP relaxations after
%applying Handelman's Theorem.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializing state space, P(x), gamma(x), nu(x), Gamma_1, and Gamma_2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

display('Initializing state space, P(x), gamma(x), nu(x), Gamma_1, and Gamma_2')
relax_type = 'LP';    
n = 2;                                                  %Dimension of state space





P.deg = 4;                                              %Degree of polynomial to be decomposed
P.coeff_total = Coeff_total(2,3);                       %Total number of coefficients of P(x)
offset_deg = 0;                                         %Degree offset

gamma.deg = P.deg+offset_deg;                           %Degree of gamma(x)
nu.deg = P.deg+offset_deg;                              %Degree of nu(x)


P.monoms = lex_exps(n,P.deg);                           %Monomials of P(x)
P.coeffs = ones(1,Coeff_total(n,P.deg));               %Set up a vector of coefficients for P(x)

gamma.monoms = lex_exps(n,gamma.deg);                   %Monomials of gamma(x)
gamma.coeff_total = Coeff_total(n,gamma.deg);            %Total number of coefficients in gamma(x)
gamma.coeffs = ones(1,Coeff_total(n,gamma.deg));        %Set up a vector of coeffs for gamma(x); these are decision variables
gamma.constraint_1 = diag(gamma.coeffs);                %Set up a constraint matrix for gamma(x)

nu.monoms = lex_exps(n,nu.deg);                         %Monomials of nu(x)
nu.coeff_total = Coeff_total(n,nu.deg);                 %Total number of coefficiens in nu(x) 
nu.coeffs = ones(1,Coeff_total(n,nu.deg));              %Set up a vector of coeffs for nu(x); these are decision variables
nu.constraint_1 = diag(nu.coeffs);                      %Set up a constraint matrix for gamma(x)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating Hessian matrices for gamma(x) and nu(x).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating Hessian matrices for gamma(x) and nu(x)...')

gammma.grad_monoms = zeros(2,Coeff_total(n,gamma.deg-1));       %Set up monomials for gradient of gamma(x).
nu.grad_monoms = zeros(2,Coeff_total(n,nu.deg-1))';             %Set up  monomials for gradient of nu(x).

gammma.hess_monoms = zeros(2,Coeff_total(n,gamma.deg-2));       %Set up monomials for gradient of gamma(x).
nu.hess_monoms = zeros(2,Coeff_total(n,nu.deg-2))';             %Set up monomials for gradient of nu(x).

%Set up gradient of gamma(x)
for i=1:n
    temp = gamma.monoms(:,i) - ones(Coeff_total(n,gamma.deg),1);
    temp_coeffs = gamma.monoms(:,i).*gamma.coeffs(:);
    temp_constr = gamma.constraint_1*diag(gamma.monoms(:,i));
    gamma.grad_monoms(:,i) = temp(temp>=0);
    gamma.grad_coeffs(:,i) = temp_coeffs(temp_coeffs>0);
    gamma.grad_constraint(:,:,i) = temp_constr(any(temp_constr,2),:);    
end


%Set up gradient of nu(x)
for i=1:n
    temp = nu.monoms(:,i) - ones(Coeff_total(n,nu.deg),1);
    temp_coeffs = nu.monoms(:,i).*nu.coeffs(:);
    temp_constr = nu.constraint_1*diag(nu.monoms(:,i));
    nu.grad_monoms(:,i) = temp(temp>=0);
    nu.grad_coeffs(:,i) = temp_coeffs(temp_coeffs>0);
    nu.grad_constraint(:,:,i) = temp_constr(any(temp_constr,2),:);  
end


%Set up Hessian of gamma(x)

for i=1:n
   for j=1:n
        temp = gamma.grad_monoms(:,i) - ones(Coeff_total(n,gamma.deg-1),1);
        temp_coeffs = gamma.grad_monoms(:,i).*gamma.grad_coeffs(:,j);
        s = gamma.grad_monoms(:,i);
        for u=1:length(s)
            temp_constraint(u,:) = gamma.grad_constraint(u,:,j)*s(u);
        end
        t = temp_constraint;
        tt=t(any(t,2),:);
        gamma.hess_monoms(:,i) = temp(temp>=0);
        gamma.hess_coeffs(i,j,:) = temp_coeffs(temp_coeffs>0);
        gamma.hess_constraint(:,:,i,j) = tt;
         
   end   
end
%Set up Hessian of nu(x)
for i=1:n
   for j=1:n
        temp = nu.grad_monoms(:,i) - ones(Coeff_total(n,nu.deg-1),1);
        temp_coeffs = nu.grad_monoms(:,i).*nu.grad_coeffs(:,j);
                s = nu.grad_monoms(:,i);
        for u=1:length(s)
            temp_constraint(u,:) = nu.grad_constraint(u,:,j)*s(u);
        end
        t = temp_constraint;
        tt=t(any(t,2),:);
        nu.hess_monoms(:,i) = temp(temp>=0);
        nu.hess_coeffs(i,j,:) = temp_coeffs(temp_coeffs>0);
        nu.hess_constraint(:,:,i,j) = tt;
         
   end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating Hessian matrices for gamma(x) and nu(x).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating quadratic forms of Hessian matrices of gamma and nu...')

%Set up y^T*H_gamma*y
gamma.quadhess_monoms = lex_exps(2*n,gamma.deg);
gamma.quadhess_coeffs = zeros(1,Coeff_total(2*n,gamma.deg)).';
gamma.quadhess_constraint = zeros(Coeff_total(2*n,gamma.deg),Coeff_total(n,gamma.deg));
for i=1:n
    for j=1:n
        for k=1:Coeff_total(n,gamma.deg-2)
            temp_exp_1 = zeros(1,n);
            temp_exp_1(i) = temp_exp_1(i)+1;
            temp_exp_1(j) = temp_exp_1(j)+1;
            gamma.hess_monoms(k,:);
            temp_exp = [temp_exp_1,gamma.hess_monoms(k,:)];
            temp_lex = lex_index_nh(temp_exp);
            temp_coeff = gamma.hess_coeffs(i,j,k);
            temp_constraint = gamma.hess_constraint(k,:,i,j);
            gamma.quadhess_coeffs(temp_lex) = gamma.quadhess_coeffs(temp_lex)+temp_coeff;
            gamma.quadhess_constraint(temp_lex,:) = temp_constraint + gamma.quadhess_constraint(temp_lex,:);
        end
    end
end


%Set up y^T*H_nu*y
nu.quadhess_monoms = lex_exps(2*n,nu.deg);
nu.quadhess_coeffs = zeros(1,Coeff_total(2*n,gamma.deg)).';
nu.quadhess_constraint = zeros(Coeff_total(2*n,nu.deg),Coeff_total(n,nu.deg));

for i=1:n
    for j=1:n        
        for k=1:Coeff_total(n,nu.deg-2)
            temp_exp_1 = zeros(1,n);
            temp_exp_1(i) = temp_exp_1(i)+1;
            temp_exp_1(j) = temp_exp_1(j)+1;
            temp_exp = [temp_exp_1,nu.hess_monoms(k,:)];
            temp_lex = lex_index_nh(temp_exp);
            temp_coeff_1 = nu.hess_coeffs(i,j,k);
            temp_constraint = nu.hess_constraint(k,:,i,j);
            nu.quadhess_coeffs(temp_lex) = nu.quadhess_coeffs(temp_lex)+ temp_coeff_1;   
            nu.quadhess_constraint(temp_lex,:) = temp_constraint + nu.quadhess_constraint(temp_lex,:);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding Handelman representations for the quadratic forms of gamma(x) and
%nu(x).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if((strcmp(relax_type,'LP')==1))
disp('Finding Handelman representations for the quadratic forms of gamma(x) and nu(x)...')


%Set up first polytope
Gamma_1 = [0 1 0 0 0; 1 -1 0 0 0; 0 0 1 0 0 ; 1 0 -1 0 0; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
n_gamma_1 = 7;
%Set up second polytope
Gamma_2 = [0 1 0 0 0; 1 -1 0 0 0; 0 0 1 0 0 ; 1 0 -1 0 0; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
n_gamma_2 = 7;

% Set up first subpolytope (gamma(x))
% Gamma_1{1} = [0 1 0 0 0; 1 -1 -1 0 0; 0 0 1 0 0 ; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
% n_gamma_1{1} = 6;
% Set up first subpolytope (nu(x))
% Gamma_2{1} = [0 1 0 0 0; 1 -1 -1 0 0; 0 0 1 0 0 ; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
% n_gamma_2{1} = 6;
% 
% Set up second subpolytope (gamma(x))
% Gamma_1{1} = [0 -1 0 0 0; 1 1 -1 0 0; 0 0 1 0 0 ; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
% n_gamma_1{1} = 6;
% Set up second subpolytope (nu(x))
% Gamma_2{1} = [0 -1 0 0 0; 1 1 -1 0 0; 0 0 1 0 0 ; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
% n_gamma_2{1} = 6;
% 
% Set up third subpolytope (gamma(x))
% Gamma_1{1} = [0 -1 0 0 0; 1 1 1 0 0; 0 0 -1 0 0 ; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
% n_gamma_1{1} = 6;
% Set up third subpolytope (nu(x))
% Gamma_2{1} = [0 -1 0 0 0; 1 1 1 0 0; 0 0 -1 0 0 ; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
% n_gamma_2{1} = 6;
% 
% Set up fourth subpolytope (gamma(x))
% Gamma_1{1} = [0 1 0 0 0; 1 -1 1 0 0; 0 0 -1 0 0 ; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
% n_gamma_1{1} = 6;
% Set up fourth subpolytope (nu(x))
% Gamma_2{1} = [0 1 0 0 0; 1 -1 1 0 0; 0 0 -1 0 0 ; 0 0 0 1 0;0 0 0 0 1; 1 0 0 -1 -1];
% n_gamma_2{1} = 6;

%Set up Handelman form for Hessian of gamma(x)
disp('Set up Handelman form for Hessian of gamma(x)... ')
n_states = 2*n;
V_gamma.monom_exps = lex_exps(2*n,gamma.deg);
V_gamma.poly_exps = lex_exps(n_gamma_1,gamma.deg);
V_gamma.poly_coeffs = Coeff_total(n_gamma_1,gamma.deg);
V_gamma.monom_coeffs = Coeff_total(2*n,gamma.deg);
V_gamma.temp_coeffs = zeros(V_gamma.monom_coeffs,V_gamma.poly_coeffs);

for i=1:V_gamma.poly_coeffs
    temp_cell = cell(n_gamma_1,1);
    tempexp = V_gamma.poly_exps(i,:);
    for j=1:n_gamma_1
        temp_cell{j} = multinomial(Gamma_1(j,:),tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    V_gamma.temp_coeffs(1:K,i) = temp_coeff;
end

%Set up Handelman form for Hessian of nu(x)
disp('Set up Handelman form for Hessian of nu(x)... ')

V_nu.monom_exps = lex_exps(2*n,nu.deg);
V_nu.poly_exps = lex_exps(n_gamma_2,nu.deg);
V_nu.poly_coeffs = Coeff_total(n_gamma_2,nu.deg);
V_nu.monom_coeffs = Coeff_total(2*n,nu.deg);
V_nu.temp_coeffs = zeros(V_nu.monom_coeffs,V_nu.poly_coeffs);

for i=1:V_nu.poly_coeffs
    temp_cell = cell(n_gamma_2,1);
    tempexp = V_nu.poly_exps(i,:);
    for j=1:n_gamma_2
        temp_cell{j} = multinomial(Gamma_2(j,:),tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    V_nu.temp_coeffs(1:K,i) = temp_coeff;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DC Decomposition using Semidefinite Programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if((strcmp(relax_type,'SDP')==1))
    disp('Proceeding to set up and solve Handelman based Semidefinite Program...')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No correct solver detected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if((strcmp(relax_type,'LP')==0)&(strcmp(relax_type,'SDP')==0))
%     
%     disp('No correct solver chosen. Please select either LP or SDP. The program will terminate now.')
%     clc
%     clear all
%     close all
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DC Decomposition using Linear Programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if((strcmp(relax_type,'LP')==1))
    disp('Proceeding to set up and solve Handelman based Linear Program...')
end

%Setting up equality cosnstraint matrix A_eq:

A_eq = [gamma.constraint_1,-nu.constraint_1,zeros(gamma.coeff_total,V_gamma.poly_coeffs),zeros(nu.coeff_total,V_nu.poly_coeffs);...
       gamma.quadhess_constraint, zeros(V_nu.monom_coeffs,nu.coeff_total), -V_gamma.temp_coeffs, zeros(V_nu.monom_coeffs, V_nu.poly_coeffs);...
       zeros(V_gamma.monom_coeffs,gamma.coeff_total),nu.quadhess_constraint, zeros(V_gamma.monom_coeffs, V_gamma.poly_coeffs),-V_nu.temp_coeffs];
   
 %Setting up equality constraint vector b_eq:

  b_eq= [zeros(nu.coeff_total - P.coeff_total,1);P.coeffs';zeros(V_nu.monom_coeffs+V_gamma.monom_coeffs,1)];
 

 %Setting up inequality constraint matrix A_ineq:
 
 A_ineq = blkdiag(zeros(gamma.coeff_total),zeros(nu.coeff_total),-eye(V_gamma.poly_coeffs),-eye(V_nu.poly_coeffs));
 
 %Setting up inequality constraint vector b_ineq:
b_ineq = zeros(gamma.coeff_total+nu.coeff_total+V_gamma.poly_coeffs+V_nu.poly_coeffs,1);

%Setting up cost function:
cost_fun = b_ineq;                          %This may be changed later if needed...
tic
%Solve the linear program with linprog:
[XX,fval,exitflag] = linprog(cost_fun, A_ineq, b_ineq, A_eq, b_eq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc



