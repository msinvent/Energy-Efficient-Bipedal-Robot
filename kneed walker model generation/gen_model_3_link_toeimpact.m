
    clear all;
    clc;
    close all;
	%% variable declarations
	%
	syms q1 q2 q3 z1 z2 real
	syms dq1 dq2 dq3 dz1 dz2 real
	syms m1 m2 Mh g L1 L2 a1 a2 real
	
	%% generalized coordinates
	%
	q=[q1;q2;q3];
	qe=[q;z1;z2];

	%% first derivative of generalized coordinates
	%
	dq=[dq1;dq2;dq3];
	dqe=[dq;dz1;dz2];

	%% Generate De matrix.  To do so, we must go through the trouble of
	%% defining the positions of the masses using the augmented
	%% configuration variables.
	%

	%% position of masses in system
	%
	p_m1=[z1+(L1-a1)*cos(q1); z2+(L1-a1)*sin(q1)];
    p_m2=[z1+(L1+L2-a2)*cos(q1) ; z2+(L1+L2-a2)*sin(q1)];
    p_Mh=[z1+(L1+L2)*cos(q1) ; z2+(L1+L2)*sin(q1)];
	p_m3=p_Mh + [a2*cos(q1+q2); a2*sin(q1+q2)];
	p_m4=p_Mh + [L2*cos(q1+q2); L2*sin(q1+q2)] + [a1*cos(q1+q2+q3); a1*sin(q1+q2+q3)];

	%% velocities of masses in system
	%
	v_m1=jacobian(p_m1,qe)*dqe;
    v_m2=jacobian(p_m2,qe)*dqe;
	v_Mh=jacobian(p_Mh,qe)*dqe;
	v_m3=jacobian(p_m3,qe)*dqe;
    v_m4=jacobian(p_m4,qe)*dqe;

	%% kinetic energy of masses in system
	%
	KE_m1=simple(m1/2*(v_m1'*v_m1));
    KE_m2=simple(m2/2*(v_m2'*v_m2));
	KE_Mh=simple(Mh/2*(v_Mh'*v_Mh));
	KE_m3=simple(m2/2*(v_m3'*v_m3));
    KE_m4=simple(m1/2*(v_m4'*v_m4));
	

	%% total kinetic energy of system
	%
	KE=KE_m1+KE_m2+KE_Mh+KE_m3+KE_m4;

	%% potential energy of masses in system
	%
	PE_m1=p_m1(2)*m1*g;
    PE_m2=p_m2(2)*m2*g;
	PE_Mh=p_Mh(2)*Mh*g;
	PE_m3=p_m3(2)*m2*g;
    PE_m4=p_m4(2)*m1*g;
	

	%% total potential energy of system
	%
	PE=PE_m1+PE_m2+PE_Mh+PE_m3+PE_m4;

	De=simplify(jacobian(jacobian(KE,dqe).',dqe))

	N=max(size(qe));
    syms Ce
	for k=1:N,
		for j=1:N,
			Ce(k,j)=0*g;
			for i=1:N,
				Ce(k,j)=Ce(k,j)+1/2*(diff(De(k,j),qe(i))+...
					diff(De(k,i),qe(j))-...
					diff(De(i,j),qe(k)))*dqe(i);
			end
		end
	end

	Ge=jacobian(PE,qe).';
    
    Psi=p_Mh + [L2*cos(q1+q2); L2*sin(q1+q2)] + [L1*cos(q1+q2+q3); L1*sin(q1+q2+q3)]
	%% Used to calculate post-impact conditions
	%
	E=jacobian(Psi,qe)


	%% Now that We're done generating the De matrix we now can use the
	%% non-augmented configuration variables to derive the required
	%% matrices.
	%

	%% position of masses in system
	%
	p_m1=[(L1-a1)*cos(q1); (L1-a1)*sin(q1)];
    p_m2=[(L1+L2-a2)*cos(q1) ; (L1+L2-a2)*sin(q1)];
    p_Mh=[(L1+L2)*cos(q1) ; (L1+L2)*sin(q1)];
	p_m3=p_Mh + [a2*cos(q1+q2); a2*sin(q1+q2)];
	p_m4=p_Mh + [L2*cos(q1+q2); L2*sin(q1+q2)] + [a1*cos(q1+q2+q3); a1*sin(q1+q2+q3)];


	%% velocities of masses in system
	%
    v_m1=jacobian(p_m1,q)*dq;
    v_m2=jacobian(p_m2,q)*dq;
	v_Mh=jacobian(p_Mh,q)*dq;
	v_m3=jacobian(p_m3,q)*dq;
    v_m4=jacobian(p_m4,q)*dq;

	%% kinetic energy of masses in system
	%
	KE_m1=simple(m1/2*(v_m1'*v_m1));
    KE_m2=simple(m2/2*(v_m2'*v_m2));
	KE_Mh=simple(Mh/2*(v_Mh'*v_Mh));
	KE_m3=simple(m2/2*(v_m3'*v_m3));
    KE_m4=simple(m1/2*(v_m4'*v_m4));
	

	%% total kinetic energy of system
	%
	KE=KE_m1+KE_m2+KE_Mh+KE_m3+KE_m4;

	%% potential energy of masses in system
	%
	PE_m1=p_m1(2)*m1*g;
    PE_m2=p_m2(2)*m2*g;
	PE_Mh=p_Mh(2)*Mh*g;
	PE_m3=p_m3(2)*m2*g;
    PE_m4=p_m4(2)*m1*g;
	

	%% total potential energy of system
	%
	PE=PE_m1+PE_m2+PE_Mh+PE_m3+PE_m4;
	
	%% the Lagragian
	%
	%Lag=KE-PE;

	%% Form D, C, G, B, and F matrices of
	%%
	%% D(q)ddq+C(q,dq)dq+G(q)=B*tau
	%%
	%% where tau=[tau1; tau2]
	%
	D=simplify(jacobian(jacobian(KE,dq).',dq))

	N=max(size(q));
    syms C
	for k=1:N,
		for j=1:N,
			C(k,j)=0*g;
			for i=1:N,
				C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+...
					diff(D(k,i),q(j))-...
					diff(D(i,j),q(k)))*dq(i);
			end
		end
    end
    C = simplify(C)
	G=jacobian(PE,q).';
    G = simplify(G)
