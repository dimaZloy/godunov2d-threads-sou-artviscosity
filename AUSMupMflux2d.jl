

# AUSM+M  flux splitting
# based on Chen et al., An improved AUSM-family scheme with robustness and accuracy for all Mach number flows 
# Applied Mathematical Modelling 77 (2020) 1065-1081


#function get_gas_epsilon(p::Float64, rho::Float64, gamma::Float64)::Float64
# return p/rho/(gamma-1.0); 
#end

# @inline  function P_m(M::Float64,AUSM_ALFA::Float64)::Float64
	# if (abs(M)>=1.0)
		# return 0.5*(1.0-sign(M));
	# else
		# return Palfa_m(M,AUSM_ALFA);
	# end
# end

# @inline  function P_p(M::Float64,AUSM_ALFA::Float64)::Float64
	# if (abs(M)>=1.0)
		# return 0.5*(1.0+sign(M));
	# else
		# return Palfa_p(M,AUSM_ALFA);
	# end
# end

# @inline  function Palfa_m(M::Float64,AUSM_ALFA::Float64)::Float64
	# return  0.25*(M-1.0)*(M-1.0)*(2.0+M)-AUSM_ALFA*M*(M*M-1.0)*(M*M-1.0);
# end

# @inline  function Palfa_p(M::Float64,AUSM_ALFA::Float64)::Float64
	# return  0.25*(M+1.0)*(M+1.0)*(2.0-M)+AUSM_ALFA*M*(M*M-1.0)*(M*M-1.0);
# end

# @inline  function Mbetta_p(M::Float64,AUSM_BETTA::Float64)::Float64
	# return  0.25*(M+1.0)*(M+1.0)+AUSM_BETTA*(M*M-1.0)*(M*M-1.0);
# end

# @inline  function Mbetta_m(M::Float64,AUSM_BETTA::Float64)::Float64
	# return -0.25*(M-1.0)*(M-1.0)-AUSM_BETTA*(M*M-1.0)*(M*M-1.0);
# end

# @inline  function M_p(M::Float64,AUSM_BETTA::Float64)::Float64
	# if (abs(M)>=1.0)
		# return 0.5*(M+abs(M));
	# else
		# return Mbetta_p(M,AUSM_BETTA);
	# end	
# end

# @inline  function M_m(M::Float64,AUSM_BETTA::Float64)::Float64
	# if (abs(M)>=1.0)
		# return  0.5*(M-abs(M));
	# else
		# return Mbetta_m(M,AUSM_BETTA);
	# end
# end



@inline  function P_m(M::Float64,AUSM_ALFA::Float64)::Float64	
	return (abs(M)>=1.0) ? 0.5*(1.0-sign(M)) : Palfa_m(M,AUSM_ALFA)
end

@inline  function P_p(M::Float64,AUSM_ALFA::Float64)::Float64
	return (abs(M)>=1.0) ? 0.5*(1.0+sign(M)) : Palfa_p(M,AUSM_ALFA)
end

@inline  function Palfa_m(M::Float64,AUSM_ALFA::Float64)::Float64
	return  0.25*(M-1.0)*(M-1.0)*(2.0+M)-AUSM_ALFA*M*(M*M-1.0)*(M*M-1.0);
end

@inline  function Palfa_p(M::Float64,AUSM_ALFA::Float64)::Float64
	return  0.25*(M+1.0)*(M+1.0)*(2.0-M)+AUSM_ALFA*M*(M*M-1.0)*(M*M-1.0);
end

@inline  function Mbetta_p(M::Float64,AUSM_BETTA::Float64)::Float64
	return  0.25*(M+1.0)*(M+1.0)+AUSM_BETTA*(M*M-1.0)*(M*M-1.0);
end

@inline  function Mbetta_m(M::Float64,AUSM_BETTA::Float64)::Float64
	return -0.25*(M-1.0)*(M-1.0)-AUSM_BETTA*(M*M-1.0)*(M*M-1.0);
end

@inline  function M_p(M::Float64,AUSM_BETTA::Float64)::Float64
	return (abs(M)>=1.0) ? 0.5*(M+abs(M)) : Mbetta_p(M,AUSM_BETTA)
end

@inline  function M_m(M::Float64,AUSM_BETTA::Float64)::Float64
	return (abs(M)>=1.0) ? 0.5*(M-abs(M)) : Mbetta_m(M,AUSM_BETTA)
end


@inline  function AUSMplusMFlux2d(
				rhoL::Float64,	_UL::Float64, 	_VL::Float64, 	PL::Float64,
				rhoR::Float64,	_UR::Float64, 	_VR::Float64,	PR::Float64, 
				nx::Float64,  ny::Float64, side::Float64, 
				gamma::Float64, flux::Array{Float64,1})


				
	VL_tilda::Float64   = _UL*nx + _VL*ny;
	VR_tilda::Float64   = _UR*nx + _VR*ny;

	TLeft::Float64      = _UL*ny - _VL*nx;
	TRight::Float64     = _UR*ny - _VR*nx;

	UMAG_L2::Float64 = (_UL*_UL + _VL*_VL);
	UMAG_R2::Float64 = (_UR*_UR + _VR*_VR);

	
	# method constants 	
	AUSM_BETTA::Float64 = 1.0/8.0;
	AUSM_ALFA::Float64  = 3.0/16.0;
	
	#htL =  get_gas_epsilon(PL,rhoL,gamma)+ 0.5*(UMAG_L2) +  PL/rhoL; 
	#htR =  get_gas_epsilon(PR,rhoR,gamma)+ 0.5*(UMAG_R2) +  PR/rhoR;

	htL::Float64 =  PL/rhoL/(gamma-1.0) + 0.5*(UMAG_L2) +  PL/rhoL; 
	htR::Float64 =  PR/rhoR/(gamma-1.0) + 0.5*(UMAG_R2) +  PR/rhoR;
	
	## h normal - original 
	htNN::Float64  = 0.5*(htL + htR - 0.5*(TLeft*TLeft + TRight*TRight));
	
	##htNN::Float64  = 0.5*(htL  - 0.5*UMAG_L2 + htR - 0.5*UMAG_R2) ;
	

	# original calculation of the numerical speed of sound from AUSM+ 
	# aL_dot::Float64 = sqrt(2.0*( gamma-1.0)/(gamma+1.0)*htL);
	# aL_tilda::Float64 = aL_dot*min(1.0,aL_dot/abs(VL_tilda));
	# aR_dot::Float64 = sqrt(2.0*( gamma-1.0)/(gamma+1.0)*htR);
	# aR_tilda::Float64 = aR_dot*min(1.0,aR_dot/abs(VR_tilda));
	# a12::Float64 = min(aL_tilda,aR_tilda);
	
	
	# compute the numerical speed of sound 
	cs::Float64 = sqrt(2.0*(gamma-1.0)/(gamma+1.0)*htNN);
	c12::Float64 = (0.5*(VL_tilda + VR_tilda) >= 0.0) ? cs*cs/max(abs(VL_tilda),cs) : cs*cs/max(abs(VR_tilda),cs) ;

	# right/left Mach numbers 
	MR::Float64 = VR_tilda/c12;
	ML::Float64 = VL_tilda/c12;
	
	rho12::Float64 = 0.5*(rhoL+rhoR);
	M::Float64 = min(1.0, max(abs(ML),abs(MR)));
	f::Float64 = 0.5*(1-cos(pi*M));
	h::Float64 = min(PL/PR, PR/PL); ### TBD 
	g::Float64 = (1.0 + cos(pi*h))*0.5;
	
	# pressure diffusive term for mass flux 
	MP::Float64 = -0.5*(1.0-f)*(PR-PL)/rho12/c12/c12*(1-g);
	

	
	KSI_P_L::Float64 = (abs(ML) < 1.0) ?  0.25*(ML+1.0)*(ML+1.0)*(2.0-ML)+AUSM_ALFA*ML*(ML*ML-1.0)*(ML*ML-1.0) : 0.5*(1.0 + sign(ML)) ;
	KSI_M_R::Float64 = (abs(MR) < 1.0) ?  0.25*(MR-1.0)*(MR-1.0)*(2.0+MR)+AUSM_ALFA*MR*(MR*MR-1.0)*(MR*MR-1.0) : 0.5*(1.0 - sign(MR)) ;
	
	
	# velocity diffusion term for pressure flux 
	
	Pux::Float64 = -g*gamma*(PL+PR)/2/c12*KSI_P_L*KSI_M_R*(_UR-_UL);
	Puy::Float64 = -g*gamma*(PL+PR)/2/c12*KSI_P_L*KSI_M_R*(_VR-_UL);
	
	# original from AUSM+  (Mp == 0)
	#m_dot12::Float64 = M_p(MLeft,AUSM_BETTA)+M_m(MRight,AUSM_BETTA);
	m_dot12::Float64 = M_p(ML,AUSM_BETTA)+M_m(MR,AUSM_BETTA) + MP ;
	
	# interface pressure 
	p12::Float64 = P_p(ML,AUSM_ALFA)*PL + P_m(MR,AUSM_ALFA)*PR;


	m_dot12_p::Float64 = 0.5*(m_dot12+abs(m_dot12));
	m_dot12_m::Float64 = 0.5*(m_dot12-abs(m_dot12));

	
	p = zeros(Float64,4);
	F_LEFT = zeros(Float64,4);
	F_RIGHT = zeros(Float64,4);
	

	#p[1] = 0.0;
	p[2] = p12*nx + Pux;
	p[3] = p12*ny + Puy; 
	#p[4] = 0.0;	
	
	
	F_LEFT[1] = rhoL;
	F_LEFT[2] = rhoL*_UL;
	F_LEFT[3] = rhoL*_VL;
	F_LEFT[4] = rhoL*htL;

	F_RIGHT[1] = rhoR;
	F_RIGHT[2] = rhoR*_UR;
	F_RIGHT[3] = rhoR*_VR;
	F_RIGHT[4] = rhoR*htR; 
	
	
	# flux[1] = -( a12*(m_dot12_p*F_LEFT[1] + m_dot12_m*F_RIGHT[1]) + p[1] )*side;
	# flux[2] = -( a12*(m_dot12_p*F_LEFT[2] + m_dot12_m*F_RIGHT[2]) + p[2] )*side;
	# flux[3] = -( a12*(m_dot12_p*F_LEFT[3] + m_dot12_m*F_RIGHT[3]) + p[3] )*side;
	# flux[4] = -( a12*(m_dot12_p*F_LEFT[4] + m_dot12_m*F_RIGHT[4]) + p[4] )*side;
	
	flux[1] = -( c12*(m_dot12_p*F_LEFT[1] + m_dot12_m*F_RIGHT[1]) + p[1] )*side;
	flux[2] = -( c12*(m_dot12_p*F_LEFT[2] + m_dot12_m*F_RIGHT[2]) + p[2] )*side;
	flux[3] = -( c12*(m_dot12_p*F_LEFT[3] + m_dot12_m*F_RIGHT[3]) + p[3] )*side;
	flux[4] = -( c12*(m_dot12_p*F_LEFT[4] + m_dot12_m*F_RIGHT[4]) + p[4] )*side;
	

	
end
