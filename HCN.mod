TITLE Hyperpolarisation activated mixed cation current (Ih)

COMMENT
Adapted from FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron
ENDCOMMENT

NEURON {
	SUFFIX HCN
	NONSPECIFIC_CURRENT i
	RANGE ghbar, eh
	GLOBAL ninf, ntau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	v	 	(mV)
	
	ghbar = .0001	(S/cm2)

	eh = -30	(mV)
 
Q10 = 3 (1) 
  Q10TEMP = 22 (degC) 
hcnfacm=1
vhhcnm=-90.1
}

ASSIGNED {
	i (mA/cm2)
	ninf
	ntau
 celsius (degC) 
  qt (1) 

}

STATE {
	n
}

INITIAL {
	rates(v)
	n = ninf
 qt = Q10^((celsius-Q10TEMP)/10) 
: qt = 1 
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = ghbar*n*(v - eh)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/ntau
}

PROCEDURE rates(v (mV)) {
	ninf = 1/(1+exp((v+vhhcnm)/9.9))
	ntau =hcnfacm*((1000 * (.19 + .72*exp(-((v-(-81.5))/11.9)^2))) / qt)
}