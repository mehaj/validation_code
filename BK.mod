TITLE BK calcium-activated potassium current

:Adapted from FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron


UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX BK
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gkbar,gk,zinf,ik
}


PARAMETER {
	celsius=37	(degC)
	v		(mV)
	gkbar=.08	(mho/cm2)	: Maximum Permeability
	cai = .04e-3	(mM)
:	ek  = -85	(mV)
	dt		(ms)
	bkfacm=1
	vhbkm=-27.91
}


ASSIGNED {
	ik		(mA/cm2)
	minf
	mexp
	zinf
	zexp
	gk
        ek (mV)
}

STATE {	m z }		: fraction of open channels

BREAKPOINT {
	SOLVE state
:	gk = gkbar*m*z*z
	ik = gkbar*m*z*z*(v - ek)
}
:UNITSOFF
:LOCAL fac

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE state() {	: exact when v held constant; integrates over dt step
	rate(v, cai)
	m = m + mexp*(minf - m)
	z = z + zexp*(zinf - z)
	VERBATIM
	return 0;
	ENDVERBATIM
}

INITIAL {
	rate(v, cai)
	m = minf
	z = zinf
}

FUNCTION alp(v (mV), ca (mM)) (1/ms) { :callable from hoc
	alp = 400/(ca*1000)
}

FUNCTION bet(v (mV)) (1/ms) { :callable from hoc
	bet = 0.11/exp((v-35)/14.9)
}

PROCEDURE rate(v (mV), ca (mM)) { :callable from hoc
	LOCAL a,b
	a = alp(v,ca)
	zinf = 1/(1+a)
	zexp = (1 - exp(-dt/10))
	b = bet(v)
	:minf = 7.5/(7.5+b)
	minf=1/(1+exp((vhbkm-v)/14.9))
	mexp = (1 - exp(-dt*(7.5+b)/bkfacm))
}
:UNITSON
