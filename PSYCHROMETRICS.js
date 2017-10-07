var err=1E-8;
var WLimit=6.6E-6;
function Altitude2AtmosphericPressure(Z){
	var P;
	if(typeof(Z)!="number") return -99999;
	P=101.325*Math.pow((1-2.25577*0.00001*Z),5.2559);
	return P;
}

function T2P_n100to0C(t){	
	var C1=-5.6745359E+3;
	var C2=6.3925247;
	var C3=-9.6778430E-3;
	var C4=6.2215701E-7;
	var C5=2.0747825E-9;
	var C6=-9.4840240E-13;
	var C7=4.1635019;	
	var p;
	var T;
	
	if(typeof(t)!="number") return -99999;
	T=t+273.15;
	p=C1/T+C2+C3*T+C4*T*T+C5*T*T*T+C6*T*T*T*T+C7*Math.log(T);
	p=Math.exp(p)/1000; //p is saturated pressure of t,kPa
	return p;
}

function T2P_0to200C(t){	
	var C8=-5.8002206E+3;
	var C9=1.3914993;
	var C10=-4.8640239E-2;
	var C11=4.1764768E-5;
	var C12=-1.4452093E-8;
	var C13=6.5459673;	
	var p;
	var T;
	
	if(typeof(t)!="number") return -99999;
	T=t+273.15;
	p=C8/T+C9+C10*T+C11*T*T+C12*T*T*T+C13*Math.log(T);
	p=Math.exp(p)/1000; //p is saturated pressure of t,kPa
	return p;
}

function T2P_Reg4(t){ //t,℃
	var n1=0.11670521452767E+4
	var n2=-0.72421316703206E+6
	var n3=-0.17073846940092E+2
	var n4=0.12020824702470E+5
	var n5=-0.32325550322333E+7
	var n6=0.14915108613530E+2
	var n7=-0.48232657361591E+4
	var n8=0.40511340542057E+6
	var n9=-0.23855557567849
	var n10=0.65017534844798E+3
	var p;
	var T;
	var theta;
	var A,B,C;
	if(typeof(t)!="number") return -99999;
	T=t+273.15;
	theta=T+n9/(T-n10);
	A=theta*theta+n1*theta+n2;
	B=n3*theta*theta+n4*theta+n5;
	C=n6*theta*theta+n7*theta+n8;
	p=Math.pow(2*C/(-B+Math.sqrt(B*B-4*A*C)),4);
	p=p*1000;//p is saturated pressure of t,kPa
	return p;	
}

function P2T_Reg4(p){  //p,kPa
	var n1=0.11670521452767E+4
	var n2=-0.72421316703206E+6
	var n3=-0.17073846940092E+2
	var n4=0.12020824702470E+5
	var n5=-0.32325550322333E+7
	var n6=0.14915108613530E+2
	var n7=-0.48232657361591E+4
	var n8=0.40511340542057E+6
	var n9=-0.23855557567849
	var n10=0.65017534844798E+3
	var t;
	var T;
	var beta;
	var D,E,F,G;
	if(typeof(p)!="number") return -99999;
	beta=Math.pow(p/1000,1/4);
	E=beta*beta+n3*beta+n6;
	F=n1*beta*beta+n4*beta+n7;
	G=n2*beta*beta+n5*beta+n8;
	D=2*G/(-F-Math.sqrt(F*F-4*E*G));
	T=(n10+D-Math.sqrt((n10+D)*(n10+D)-4*(n9+n10*D)))/2;
	t=T-273.15;//t is saturated temperature of p,℃
	return t;	
}

function T2P(t){ //t,℃
	var p;
	if(typeof(t)!="number") return -99999;
	if(t>=-100 && t<0)
		return T2P_n100to0C(t);
	else if(t>=0 && t<=200)
		return T2P_0to200C(t);
	else if(t>200)
		return T2P_Reg4(t);
	else
		return -99999;
}

function P2T(p){//p,kPa
	var p200=T2P(200);
	var p0=T2P(0);
	var t0,t1,t,p0,p1,p2;
	if(typeof(p)!="number") return -99999;
	if(p>p200) return P2T_Reg4(p);
	if(p>=p0){
		 //以下采用弦截法逼近求取温度
		t0 = 1;
		t1=199;
		p0=T2P_0to200C(t0);
		p1=T2P_0to200C(t1);
		t = t1-(p1-p)*(t0-t1)/(p0-p1);
        p2=T2P_0to200C(t);
		while(Math.abs(p - p2) > err){
			t0 = t1;
			t1 = t;
			p0 = p1;
			p1 = p2;
			t = t1-(p1-p)*(t0-t1)/(p0-p1);
			p2=T2P_0to200C(t);
		}
		return t;
	}
	// t:-100to0
	//以下采用弦截法逼近求取温度
	t0 = -1;
	t1=-60;
	p0=T2P_n100to0C(t0);
	p1=T2P_n100to0C(t1);
	t = t1-(p1-p)*(t0-t1)/(p0-p1);
	p2=T2P_n100to0C(t);
	while(Math.abs(p - p2) > err){
		t0 = t1;
		t1 = t;
		p0 = p1;
		p1 = p2;
		t = t1-(p1-p)*(t0-t1)/(p0-p1);
		p2=T2P_n100to0C(t);
	}
	return t;
}

function Tdew_0to93C(pw){//the dew-point temperature range of 0 to 93°C,pw is water vapor partial pressure, kPa
	var C14 = 6.54;
	var C15 = 14.526;
	var C16 = 0.7389;
	var C17 = 0.09486;
	var C18 = 0.4569;
	var alpha;
	var tdew;
	if(typeof(pw)!="number") return -99999;
	alpha=Math.log(pw);
	tdew=C14+C15*alpha+C16*alpha*alpha+C17*alpha*alpha*alpha+C18*Math.pow(pw,0.1984);
	return tdew;	
}

function Tdew_below0C(pw){//the dew-point temperature below 0°C,pw is water vapor partial pressure, kPa
	var alpha;
	var tdew;
	if(typeof(pw)!="number") return -99999;
	alpha=Math.log(pw);
	tdew=6.09+12.608*alpha+0.4959*alpha*alpha;
	return tdew;	
}

function Tdew(pw){//the dew-point temperature,pw is water vapor partial pressure, kPa
    var ts;
	if(typeof(pw)!="number") return -99999;
	if(pw>=0.6112 && pw<=T2P(93.001))
		return Tdew_0to93C(pw);
	else if(pw>T2P(93.001))
		return -99999;
	else 
		return Tdew_below0C(pw);	
}

function MoistAir_ttwp(t,tw,p){ //input:Dry-bulb temperature t,C; Wet-bulb temperature tw,C; Pressure p,kPa
	var pws_t,pws_tw,Ws_tw,W,Ws_t,mu,phi,nu,h,pw,td;
	var Rda=287.055;
	if(typeof(t)!="number" || typeof(tw)!="number" ||typeof(p)!="number") {
		var info={status:-1,t:-99999,tw:-99999,td:-99999,pw:-99999,phi:-99999,W:-99999,h:-99999,mu:-99999,nu:-99999};
		return info;
	}		
	pws_t=T2P(t);
	pws_tw=T2P(tw);
	Ws_t=0.62198*pws_t/(p-pws_t);
	Ws_tw=0.62198*pws_tw/(p-pws_tw);
	W=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);
	mu=W/Ws_t;
	phi=mu/(1-(1-mu)*(pws_t/p));
	nu=Rda*(t+273.15)*(1+1.6078*W)/p;
	h=1.006*t+W*(2501+1.805*t);
	pw=(p*W)/(0.62198+W);
	td=P2T(pw);
	W=W*1000; //W,unit:g;
	phi*=100;//phi,unit:%
	var info={status:0,t:t,tw:tw,td:td,pw:pw,phi:phi,W:W,h:h,mu:mu,nu:nu};
	return info;
}

function MoistAir_ttdp(t,td,p){//input:Dry-bulb temperature t,C; dew-point temperature td,C; Pressure p,kPa
	var pws_t,pw,W,Ws_t,mu,phi,nu,h,tw;
	var Rda=287.055;
	if(typeof(t)!="number" || typeof(td)!="number" ||typeof(p)!="number") {
		var info={status:-1,t:-99999,tw:-99999,td:-99999,pw:-99999,phi:-99999,W:-99999,h:-99999,mu:-99999,nu:-99999};
		return info;
	}		
	pws_t=T2P(t);
	pw=T2P(td);//partial pressure,or Saturation pressure of td
	Ws_t=0.62198*pws_t/(p-pws_t);
	W=0.62198*pw/(p-pw);
	mu=W/Ws_t;
	phi=mu/(1-(1-mu)*(pws_t/p));
	nu=Rda*(t+273.15)*(1+1.6078*W)/p;
	h=1.006*t+W*(2501+1.805*t);
	//以下采用弦截法逼近求取湿球温度tw 
	var tw0=td+0.5;
	var tw1=t-0.5;
	var pws_tw0=T2P(tw0);
	var pws_tw1=T2P(tw1);
	var Ws_tw0=0.62198*pws_tw0/(p-pws_tw0);
	var Ws_tw1=0.62198*pws_tw1/(p-pws_tw1);
	var W0=((2501-2.381*tw0)*Ws_tw0-1.006*(t-tw0))/(2501+1.805*t-4.186*tw0);
	var W1=((2501-2.381*tw1)*Ws_tw1-1.006*(t-tw1))/(2501+1.805*t-4.186*tw1);
	tw=tw1-(W1-W)*(tw0-tw1)/(W0-W1);
	var pws_tw=T2P(tw);
	var Ws_tw=0.62198*pws_tw/(p-pws_tw);
	var Wtw=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);
	while(Math.abs(W-Wtw)>err){  //需要优化
		tw0=tw1;
		W0=W1;
		tw1=tw;
		pws_tw1=T2P(tw1);
		Ws_tw1=0.62198*pws_tw1/(p-pws_tw1);
		W1=((2501-2.381*tw1)*Ws_tw1-1.006*(t-tw1))/(2501+1.805*t-4.186*tw1);
		tw=tw1-(W1-W)*(tw0-tw1)/(W0-W1);
		pws_tw=T2P(tw);
		Ws_tw=0.62198*pws_tw/(p-pws_tw);
		Wtw=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);		
	}
	//弦截法结束
	W=W*1000; //W,unit:g;
	phi*=100;//phi,unit:%
	var info={status:0,t:t,tw:tw,td:td,pw:pw,phi:phi,W:W,h:h,mu:mu,nu:nu};
	return info;
}

function MoistAir_tphip(t,phi,p){//input:Dry-bulb temperature t,C; Relative humidity,phi,%; Pressure p,kPa
	var pws_t,pw,W,Ws_t,mu,phi,nu,h,td,tw;
	var Rda=287.055;
	if(typeof(t)!="number" || typeof(phi)!="number" ||typeof(p)!="number") {
		var info={status:-1,t:-99999,tw:-99999,td:-99999,pw:-99999,phi:-99999,W:-99999,h:-99999,mu:-99999,nu:-99999};
		return info;
	}
	phi/=100;	
	pws_t=T2P(t);
	pw=pws_t*phi;//partial pressure
	Ws_t=0.62198*pws_t/(p-pws_t);
	W=0.62198*pw/(p-pw);
	mu=W/Ws_t;
	td=P2T(pw);
	nu=Rda*(t+273.15)*(1+1.6078*W)/p;
	h=1.006*t+W*(2501+1.805*t);
	//以下采用弦截法逼近求取湿球温度tw
	var tw0=td+0.5;
	var tw1=t-0.5;
	var pws_tw0=T2P(tw0);
	var pws_tw1=T2P(tw1);
	var Ws_tw0=0.62198*pws_tw0/(p-pws_tw0);
	var Ws_tw1=0.62198*pws_tw1/(p-pws_tw1);
	var W0=((2501-2.381*tw0)*Ws_tw0-1.006*(t-tw0))/(2501+1.805*t-4.186*tw0);
	var W1=((2501-2.381*tw1)*Ws_tw1-1.006*(t-tw1))/(2501+1.805*t-4.186*tw1);
	tw=tw1-(W1-W)*(tw0-tw1)/(W0-W1);
	var pws_tw=T2P(tw);
	var Ws_tw=0.62198*pws_tw/(p-pws_tw);
	var Wtw=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);
	while(Math.abs(W-Wtw)>err){
		tw0=tw1;
		W0=W1;
		tw1=tw;
		pws_tw1=T2P(tw1);
		Ws_tw1=0.62198*pws_tw1/(p-pws_tw1);
		W1=((2501-2.381*tw1)*Ws_tw1-1.006*(t-tw1))/(2501+1.805*t-4.186*tw1);
		tw=tw1-(W1-W)*(tw0-tw1)/(W0-W1);
		pws_tw=T2P(tw);
		Ws_tw=0.62198*pws_tw/(p-pws_tw);
		Wtw=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);		
	}
	//弦截法结束
	W=W*1000; //W,unit:g;
	phi*=100;//phi,unit:%
	var info={status:0,t:t,tw:tw,td:td,pw:pw,phi:phi,W:W,h:h,mu:mu,nu:nu};
	return info;
}

function TPhiP2W(t,phi,p){//input:Dry-bulb temperature t,C; Relative humidity,phi,%; Pressure p,kPa;output:W,g
	if(typeof(t)!="number" || typeof(phi)!="number" || typeof(p)!="number") return -99999;
	if(phi<0 || phi>100) return -99999;
	var pws_t,pw,W;
	phi/=100;	
	pws_t=T2P(t);
	pw=pws_t*phi;//partial pressure
	W=0.62198*pw/(p-pw);
	W*=1000;
	phi*=100;
	return W;
}

function TTwP2W(t,tw,p){//input:Dry-bulb temperature t,C; wet-bulb temperature tw,C;Pressure p,kPa;output:W,g
	if(typeof(t)!="number" || typeof(tw)!="number"|| typeof(p)!="number") return -99999;
	if(tw>t) return -99999;
	var pws_t,pws_tw,Ws_tw,W;
	pws_t=T2P(t);
	pws_tw=T2P(tw);
	Ws_tw=0.62198*pws_tw/(p-pws_tw);
	W=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);
	W*=1000;
	return W;
}

function THP2W(t,h,p){//input:Dry-bulb temperature t,C; enthalpy h,kJ/kg;Pressure p,kPa;output:W,g
	if(typeof(t)!="number" || typeof(h)!="number"|| typeof(p)!="number") return -99999;
	//if(h<0) return -99999;
	W=(h-1.006*t)/(2501+1.805*t);
	W*=1000;
	return W;	
}

function WPhiP2T(W,phi,p){//input:Humidity ratio,W,g; Relative humidity,phi,%; Pressure p,kPa;output:Dry-bulb temperature t,C
	if(typeof(W)!="number" || typeof(phi)!="number"|| typeof(p)!="number") return -99999;
	if(phi<0 || phi>100) return -99999;
	if(W<WLimit) return -99999;
	var t0,t1,t,pw0,pw1,pw,W0,W1,W2;
	phi/=100;
	W/=1000;
	//以下采用弦截法逼近求取干球温度t
	t0=5;
	t1=50;
	pw0=T2P(t0)*phi;//partial pressure
	pw1=T2P(t1)*phi;
	W0=0.62198*pw0/(p-pw0);
	W1=0.62198*pw1/(p-pw1);
	t=t1-(W1-W)*(t0-t1)/(W0-W1);
	pw=T2P(t)*phi;
	W2=0.62198*pw/(p-pw);
	while(Math.abs(W-W2)>err){
		t0=t1;
		W0=W1;
		t1=t;
		W1=W2;
		t=t1-(W1-W)*(t0-t1)/(W0-W1);
		pw=T2P(t)*phi;
		W2=0.62198*pw/(p-pw);	
	}
	phi*=100;
	W*=1000;
	return t;
}

function WTwP2T(W,tw,p){//input:Humidity ratio,W,g; wet-bulb temperature tw,C; Pressure p,kPa;output:Dry-bulb temperature t,C
	if(typeof(W)!="number" || typeof(tw)!="number"|| typeof(p)!="number") return -99999;
	if(W<WLimit) return -99999;
	var pws_tw,Ws_tw;
	var t0,t1,W0,W1,W2;
	W/=1000;	
	pws_tw=T2P(tw);
	Ws_tw=0.62198*pws_tw/(p-pws_tw);
	//以下采用弦截法逼近求取干球温度t
	t0=tw;
	t1=tw+10;
	W0=((2501-2.381*tw)*Ws_tw-1.006*(t0-tw))/(2501+1.805*t0-4.186*tw);
	W1=((2501-2.381*tw)*Ws_tw-1.006*(t1-tw))/(2501+1.805*t1-4.186*tw);
	t=t1-(W1-W)*(t0-t1)/(W0-W1);
	W2=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);
	while(Math.abs(W-W2)>err){
		t0=t1;
		W0=W1;
		t1=t;
		W1=W2;
		t=t1-(W1-W)*(t0-t1)/(W0-W1);
		W2=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);
	}
	W*=1000;
	return t;
}

function WHP2T(W,h,p){//input: Humidity ratio,W,g;enthalpy h,kJ/kg;Pressure p,kPa;output:Dry-bulb temperature t,C;
	if(typeof(W)!="number" || typeof(h)!="number"|| typeof(p)!="number") return -99999;
	//if(h<0) return -99999;
	var t;
	W/=1000;
	t=(h-2501*W)/(1.006+1.805*W);
	W*=1000;
	return t;	
}

function TWP2Phi(t,W,p){//input:Dry-bulb temperature t,C; Humidity ratio,W,g; Pressure p,kPa;output:Relative humidity,phi,%;
	if(typeof(t)!="number" || typeof(W)!="number" || typeof(p)!="number") return -99999;
	if(W<0) return -99999;
	var pws_t,pw,phi;
	W/=1000;
	pw=W*p/(0.62198+W);
	pws_t=T2P(t);
	phi=pw/pws_t;
	W*=1000;
	phi*=100;
	return phi;
}

function TWP2Tw(t,W,p){//input:Dry-bulb temperature t,C; Humidity ratio,W,g; Pressure p,kPa;output:wet-bulb temperature tw;
	if(typeof(t)!="number" || typeof(W)!="number" || typeof(p)!="number") return -99999;
	var pws_tw,Ws_tw,tw0,tw1,tw2,W0,W1,W2,tw;
	W/=1000;
	tw0=t-5;
	pws_tw=T2P(tw0);
	Ws_tw=0.62198*pws_tw/(p-pws_tw);
	W0=((2501-2.381*tw0)*Ws_tw-1.006*(t-tw0))/(2501+1.805*t-4.186*tw0);
	tw1=t-1;
	pws_tw=T2P(tw1);
	Ws_tw=0.62198*pws_tw/(p-pws_tw);
	W1=((2501-2.381*tw1)*Ws_tw-1.006*(t-tw1))/(2501+1.805*t-4.186*tw1);
	tw=tw1-(W1-W)*(tw0-tw1)/(W0-W1);
	pws_tw=T2P(tw);
	Ws_tw=0.62198*pws_tw/(p-pws_tw);
	W2=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);
	while(Math.abs(W-W2)>err){
		tw0=tw1;
		W0=W1;
		tw1=tw;
		W1=W2;
		tw=tw1-(W1-W)*(tw0-tw1)/(W0-W1);
		pws_tw=T2P(tw);
		Ws_tw=0.62198*pws_tw/(p-pws_tw);
		W2=((2501-2.381*tw)*Ws_tw-1.006*(t-tw))/(2501+1.805*t-4.186*tw);
	}	
	W*=1000;
	return tw;
}

function TWP2H(t,W,p){//input:Dry-bulb temperature t,C; Humidity ratio,W,g; Pressure p,kPa;output:enthalpy h,kJ/kg;
	if(typeof(t)!="number" || typeof(W)!="number" || typeof(p)!="number") return -99999;
	var h;
	W/=1000;
	h=1.006*t+W*(2501+1.805*t);
	W*=1000;
	return h;
}

function PhiHP2T(phi,h,p){//input: Relative humidity,phi,%; enthalpy h,kJ/kg;Pressure p,kPa;Dry-bulb temperature t,C;
	if(typeof(phi)!="number" || typeof(h)!="number" || typeof(p)!="number") return -99999;
	var moist,t,t0,t1,t2,h0,h1,h2;
	t0=10;
	t1=40;
	moist=MoistAir_tphip(t0,phi,p);
	h0=moist.h;
	moist=MoistAir_tphip(t1,phi,p);
	h1=moist.h;
	t=t1-(h1-h)*(t0-t1)/(h0-h1);
	moist=MoistAir_tphip(t,phi,p);
	h2=moist.h;
	while(Math.abs(h-h2)>err){
		t0=t1;
		t1=t;
		h0=h1;
		h1=h2;
		t=t1-(h1-h)*(t0-t1)/(h0-h1);
		moist=MoistAir_tphip(t,phi,p);
		h2=moist.h;		
	}
	return t;
}

function PhiWP2T(phi,W,p){//input: Relative humidity,phi,%; Humidity ratio,W,g;Pressure p,kPa;Dry-bulb temperature t,C;
	if(typeof(phi)!="number" || typeof(W)!="number" || typeof(p)!="number") return -99999;
	var moist,t,t0,t1,t2,W0,W1,W2;
	t0=10;
	t1=40;
	moist=MoistAir_tphip(t0,phi,p);
	W0=moist.W;
	moist=MoistAir_tphip(t1,phi,p);
	W1=moist.W;
	t=t1-(W1-W)*(t0-t1)/(W0-W1);
	moist=MoistAir_tphip(t,phi,p);
	W2=moist.W;
	while(Math.abs(W-W2)>err){
		t0=t1;
		t1=t;
		W0=W1;
		W1=W2;
		t=t1-(W1-W)*(t0-t1)/(W0-W1);
		moist=MoistAir_tphip(t,phi,p);
		W2=moist.W;		
	}
	return t;
}

function IsoPhiCurve(svg,phi,Patm,curveId,t_min,t_max,W_min,W_max,phi_min,phi_max,phi_step){
	var phi0=TWP2Phi(t_max,W_min,Patm);
	var phi1=TWP2Phi(t_min,W_min,Patm);
	var phi2=TWP2Phi(t_max,W_max,Patm);
	var t1=WPhiP2T(W_min,phi,Patm);
	if(t1<t_min){
		t1=t_min;
	}
	var t2=t_max;
	if(phi2<phi && phi<=100){
		t2=WPhiP2T(W_max,phi,Patm);
	}
	if(phi<=100){
		var t=d3.range(t1,t2-0.1,t_step);
		t.push(t2);
		var data=t.map(function(d){return {x:d,y:TPhiP2W(d,phi,Patm)};
							});						
		svg.selectAll('#'+curveId)
				.data([])
				.exit()
				.remove('path');
		svg.selectAll('#'+curveId)
				.data([data])
				.enter()
				.append("path")
				.attr('id',''+curveId)
				.attr("class","line_phi")
				.attr("d",line);
	}
return 0;
}

function IsoPhiCluster(svg,Patm,t_min,t_max,W_min,W_max,phi_min,phi_max,phi_step){	
	var phi0=TWP2Phi(t_max,W_min,Patm);
	var phi1=TWP2Phi(t_min,W_min,Patm);
	var phi2=TWP2Phi(t_max,W_max,Patm);
	var phi=d3.range(phi0,phi_max-0.01,phi_step);
	phi.push(phi_max);
	phi.forEach(function(x2){
						var t1=WPhiP2T(W_min,x2,Patm);
						if(t1<t_min){
							t1=t_min;
						}
						var t2=t_max;
						if(phi2<x2 && x2<=100){
							t2=WPhiP2T(W_max,x2,Patm);
						}
						if(x2<=100){
							var t=d3.range(t1,t2-0.1,t_step);
							t.push(t2);
							var data=t.map(function(d){return {x:d,y:TPhiP2W(d,x2,Patm)};
												});
							svg.append("path")
									.datum(data)
									.attr("class","line_phi")
									.attr("d",line);
						}
				});
return 0;
}

function IsoTwCurve(svg,tw,Patm,curveId,t_min,t_max,W_min,W_max,tw_min,tw_max,tw_step){
	var tw2=TWP2Tw(t_max,W_max,Patm);
	var tw1=WPhiP2T(W_max,100,Patm);
	var tw0=TWP2Tw(t_max,W_min,Patm);	
	var t1=tw;
	if(tw1<tw && tw<=tw2){
		t1=WTwP2T(W_max,tw,Patm);
		t2=t_max;						
	}else if(tw0<tw && tw<=tw1){
		t2=t_max;						t
	}else if(t_min<=tw && tw<tw0){
		t2=WTwP2T(W_min,tw,Patm);
	}
	if(t_min<=tw && tw<=tw2){						
		var t=d3.range(t1,t2-0.1,t_step);
		t.push(t2);
		var data=t.map(function(d){return {x:d,y:TTwP2W(d,tw,Patm)};
							});						
		svg.selectAll('#'+curveId)
				.data([])
				.exit()
				.remove('path');
		svg.selectAll('#'+curveId)
				.data([data])
				.enter()
				.append("path")
				.attr('id',''+curveId)
				.attr("class","line_tw")
				.attr("d",line);
	};
	return 0;
}

function IsoTwCluster(svg,Patm,t_min,t_max,W_min,W_max,tw_min,tw_max,tw_step){
	var tw2=TWP2Tw(t_max,W_max,Patm);
	var tw1=WPhiP2T(W_max,100,Patm);
	var tw0=TWP2Tw(t_max,W_min,Patm);						
	var tw=d3.range(tw_min,tw2+tw_step,tw_step);	
	tw.forEach(function(x2){
					var t1=x2;
					if(tw1<x2 && x2<=tw2){
						t1=WTwP2T(W_max,x2,Patm);
						t2=t_max;						
					}else if(tw0<x2 && x2<=tw1){
						t2=t_max;						t
					}else if(t_min<=x2 && x2<tw0){
						t2=WTwP2T(W_min,x2,Patm);
					}
					if(t_min<=x2 && x2<=tw2){						
						var t=d3.range(t1,t2-0.1,t_step);
						t.push(t2);
						var data=t.map(function(d){return {x:d,y:TTwP2W(d,x2,Patm)};
											});
						svg.append("path")
								.datum(data)
								.attr("class","line_tw")
								.attr("d",line);
					};
			});
	return 0;
}

function IsoHCurve(svg,h,Patm,curveId,t_min,t_max,W_min,W_max,h_min,h_max,h_step){	
	var moist=MoistAir_tphip(t_min,100,Patm);
	var h0=moist.h;
	var tt=PhiWP2T(100,W_max,Patm);
	var h1=TWP2H(tt,W_max,Patm);
	var h2=TWP2H(t_max,W_max,Patm);
	if(h1<h && h<=h2){
		var t1=WHP2T(W_max,h,Patm);
		var t2=t_max;
	}else if(h0<h & h<=h1){
		var t1=PhiHP2T(100,h,Patm);
		var t2=WHP2T(W_min,h,Patm);
		if(t2>t_max) t2=t_max;
	}else if(h<h0){
		t1=t_min;
		t2=WHP2T(W_min,h,Patm);
	}
	if(h<=h2){
		var t=d3.range(t1,t2-0.1,t_step);
		t.push(t2);
		var data=t.map(function(d){return {x:d,y:THP2W(d,h,Patm)};
						});						
		svg.selectAll('#'+curveId)
				.data([])
				.exit()
				.remove('path');
		svg.selectAll('#'+curveId)
				.data([data])
				.enter()
				.append("path")
				.attr('id',''+curveId)
				.attr("class","line_h")
				.attr("d",line);
	}
	return 0;
}


function IsoHCluster(svg,Patm,t_min,t_max,W_min,W_max,h_min,h_max,h_step){	
	var moist=MoistAir_tphip(t_min,100,Patm);
	var h0=moist.h;
	var tt=PhiWP2T(100,W_max,Patm);
	var h1=TWP2H(tt,W_max,Patm);
	var h2=TWP2H(t_max,W_max,Patm);			
	var h=d3.range(h_min,h2+h_step,h_step);
	h.forEach(function(x2){
						if(h1<x2 && x2<=h2){
							var t1=WHP2T(W_max,x2,Patm);
							var t2=t_max;
						}else if(h0<x2 & x2<=h1){
							var t1=PhiHP2T(100,x2,Patm);
							var t2=WHP2T(W_min,x2,Patm);
							if(t2>t_max) t2=t_max;
						}else if(x2<h0){
							t1=t_min;
							t2=WHP2T(W_min,x2,Patm);
						}
						if(x2<=h2){
							var t=d3.range(t1,t2-0.1,t_step);
							t.push(t2);
							var data=t.map(function(d){return {x:d,y:THP2W(d,x2,Patm)};
											});
							
							svg.append("path")
									.datum(data)
									.attr("class","line_h")
									.attr("d",line);
						}
					});
	return 0;
}


