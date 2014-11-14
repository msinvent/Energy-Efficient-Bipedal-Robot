function [] =  gen_model_knee_impact
    clear all;
    clc;
    close all;
	%% variable declarations
	%
	syms q1 q2 q3 real
	syms dq1_n dq2_n dq3_n dq1_p dq2_p
	syms m1 m2 Mh g L1 L2 a1 a2 real
	
	%% generalized coordinates
	%
	q_n=[q1;q2;q3];
	dq_n=[dq1_n;dq2_n;dq3_n];
	


	%% Now that We're done generating the De matrix we now can use the
	% non-augmented configuration variables to derive the required
	% matrices

	%% position of masses before impact
	%
	p_m1_n=[(L1-a1)*cos(q1); (L1-a1)*sin(q1)];
    p_m2_n=[(L1+L2-a2)*cos(q1) ; (L1+L2-a2)*sin(q1)];
    p_Mh_n=[(L1+L2)*cos(q1) ; (L1+L2)*sin(q1)];
	p_m3_n=p_Mh_n + [a2*cos(q1+q2); a2*sin(q1+q2)];
	p_m4_n=p_Mh_n + [L2*cos(q1+q2); L2*sin(q1+q2)] + [a1*cos(q1+q2+q3); a1*sin(q1+q2+q3)];
	% velocities of masses before 	
    v_m1_n=jacobian(p_m1_n,q_n)*dq_n;
    v_m2_n=jacobian(p_m2_n,q_n)*dq_n;
	v_Mh_n=jacobian(p_Mh_n,q_n)*dq_n;
	v_m3_n=jacobian(p_m3_n,q_n)*dq_n;
    v_m4_n=jacobian(p_m4_n,q_n)*dq_n;    
    
    
    %% position of masses after impact
	% 2 link
    q_p = [q1;q2];
    dq_p = [dq1_p;dq2_p];
	p_m1_p=[(L1-a1)*cos(q1); (L1-a1)*sin(q1)];
    p_m2_p=[(L1+L2-a2)*cos(q1) ; (L1+L2-a2)*sin(q1)];
    p_Mh_p=[(L1+L2)*cos(q1) ; (L1+L2)*sin(q1)];
	p_m3_p=p_Mh_p + [a2*cos(q1+q2); a2*sin(q1+q2)];
	p_m4_p=p_Mh_p + [(L2+a1)*cos(q1+q2); (L2+a1)*sin(q1+q2)];

	% velocities of masses after impact
	%
	v_m1_p=jacobian(p_m1_p,q_p)*dq_p;
    v_m2_p=jacobian(p_m2_p,q_p)*dq_p;
	v_Mh_p=jacobian(p_Mh_p,q_p)*dq_p;
	v_m3_p=jacobian(p_m3_p,q_p)*dq_p;
    v_m4_p=jacobian(p_m4_p,q_p)*dq_p;

%% angular momentum conservation about the hip
    sigma_3 = m2*wp((p_m3_n-p_Mh_n),v_m3_n);
    sigma_4 = m1*wp((p_m4_n-p_Mh_n),v_m4_n);
    sigma_hip_n = simplify(sigma_3 + sigma_4)
    sigma_3 = m2*wp((p_m3_p-p_Mh_p),v_m3_p);
    sigma_4 = m1*wp((p_m4_p-p_Mh_p),v_m4_p);
    sigma_hip_p = simplify(sigma_3 + sigma_4)
 %angular momentum conservation about the pivot
    sigma_1 = m1*wp((p_m3_n),v_m3_n);
    sigma_2 = m2*wp((p_m4_n),v_m4_n);
    sigma_Mh = Mh*wp((p_Mh_n),v_Mh_n);
    sigma_3 = m2*wp((p_m3_n),v_m3_n);
    sigma_4 = m1*wp((p_m4_n),v_m4_n);
    sigma_pivot_n = sigma_1 + sigma_2 + sigma_Mh + sigma_3 + sigma_4
    sigma_1 = m1*wp((p_m3_p),v_m3_p);
    sigma_2 = m2*wp((p_m4_p),v_m4_p);
    sigma_Mh = Mh*wp((p_Mh_p),v_Mh_p);
    sigma_3 = m2*wp((p_m3_p),v_m3_p);
    sigma_4 = m1*wp((p_m4_p),v_m4_p);
    sigma_pivot_p = sigma_1 + sigma_2 + sigma_Mh + sigma_3 + sigma_4
    
%% variables to solve dq1_p and dq2_p
    eqn_1 = sigma_hip_n - sigma_hip_p;
    eqn_2 = sigma_pivot_n - sigma_pivot_p;
    temp = simplify(solve(eqn_1,dq2_p));
    eqn_2 = subs(eqn_2,dq2_p,temp);
    dq1_p = simplify(solve(eqn_2,dq1_p))
    dq2_p = simplify(subs(temp,dq1_p))
    
    end

    function val = wp(x,y)
    val = x'*[0 1;-1 0]*y;
    end