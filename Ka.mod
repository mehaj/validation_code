TITLE Fast inactivating potassium current


 
COMMENT
Adapted from FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Ka
	USEION k READ ek WRITE ik
        RANGE  gkbar, gk, minf, hinf, mexp, hexp, ik
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar	= .015 (mho/cm2)
     :   ek	= -85 (mV)
	kafacm=1
	kafach=1
	vhkam=-39.71
	vhkah=-59.56
}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf hinf mexp hexp 
        ek (mV)
}
 
BREAKPOINT {
        SOLVE states
        gk = gkbar *m*m*m* m*h 
	ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

PROCEDURE states() {  :Computes state variables m, h
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE minf, mexp, hinf, hexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                :"m" potassium activation system
        alpha = 1.4/(1+exp((v+27)/(-12)))
        beta =  0.49/(1+exp((v+30)/4))
        sum = alpha + beta
        :minf = alpha/sum
	minf=1/(1+exp((vhkam-v)/9.48))
        mexp = 1 - exp(tinc*sum/kafacm)
                :"h" potassium inactivation system
        alpha = 0.0175/(1+exp((v+50)/8))
        beta = 1.3/(1+exp((v+13)/(-10)))
        sum = alpha + beta
        :hinf = alpha/sum
	hinf=1/(1+exp((v-vhkah)/7.62))

        hexp = 1 - exp(tinc*sum/kafach)
}

 
UNITSON

