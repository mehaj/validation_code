TITLE Squid axon potassium channel

 
COMMENT
 Adapted from FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Khh
        USEION k READ ek WRITE ik
        RANGE   gk,  gkbar, ik
        GLOBAL  ninf, nexp
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gkbar (mho/cm2)
     :   ek = -85(mV)
	khhfacm=1
	vhkhhm=-53
}
 
STATE {
         n
}
 
ASSIGNED {
        ik (mA/cm2)
        gk ninf nexp
         ek (mV)
}
 
BREAKPOINT {
        SOLVE states
        gk  = gkbar*n*n*n*n

        ik = gk*(v - ek)      
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	n = ninf
}

PROCEDURE states() {  :Computes state variable n 
        rates(v)      :             at the current v and dt.
        n = n + nexp*(ninf-n)
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE ninf, nexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                :"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
        sum = alpha + beta
        :ninf = alpha/sum
	ninf=1/(1+exp((vhkhhm-v)/17.14))
        nexp = 1 - exp(tinc*sum/khhfacm)
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

