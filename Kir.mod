TITLE Potassium inward rectifier current

COMMENT
adapted from
	Masoli S, Solinas S, D'Angelo E (2015) Action potential processing in a detailed Purkinje cell model reveals a critical role for axonal compartmentalization. Front Cell Neurosci 
	mammalian Purkinje neuron model

Suffix from Ubc_Kir to Kir2_3
	
ENDCOMMENT
 
NEURON { 
	SUFFIX Kir
	USEION k READ ek WRITE ik 
	RANGE gkbar, ik, g, alpha_d, beta_d, ek
	RANGE Aalpha_d, Kalpha_d, V0alpha_d
	RANGE Abeta_d, Kbeta_d, V0beta_d
	RANGE d_inf, tau_d 
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	Aalpha_d = 0.13289 (/ms)

	
	Kalpha_d = -24.3902 (mV)

	V0alpha_d = -83.94 (mV)
	Abeta_d = 0.16994 (/ms)

	
	Kbeta_d = 35.714 (mV)

	V0beta_d = -83.94 (mV)
	v (mV) 
	gkbar = 0.0009 (mho/cm2) 
	ek (mV) 
	celsius = 30 (degC) 
	kirfacm=1
	vhkirm=-65.02
} 

STATE { 
	d 
} 

ASSIGNED { 
	ik (mA/cm2) 
	d_inf 
	tau_d (ms) 
	g (mho/cm2) 
	alpha_d (/ms) 
	beta_d (/ms) 
} 
 
INITIAL { 
	rate(v) 
	d = d_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit
	g = gkbar*d   
	ik = g*(v - ek) 
	alpha_d = alp_d(v) 
	beta_d = bet_d(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	d' =(d_inf - d)/tau_d 
} 
 
FUNCTION alp_d(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	alp_d = Q10*Aalpha_d*exp((v-V0alpha_d)/Kalpha_d) 
} 
 
FUNCTION bet_d(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	bet_d = Q10*Abeta_d*exp((v-V0beta_d)/Kbeta_d) 
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_d, b_d 
	TABLE d_inf, tau_d  
	DEPEND Aalpha_d, Kalpha_d, V0alpha_d, 
	       Abeta_d, Kbeta_d, V0beta_d, celsius FROM -100 TO 100 WITH 200 
	a_d = alp_d(v)  
	b_d = bet_d(v) 
	tau_d = kirfacm*(1/(a_d + b_d)) 
	:d_inf = a_d/(a_d + b_d) 
	d_inf=1/(1+exp((vhkirm-v)/76.92))
} 

