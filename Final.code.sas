************************************************************;
************************************************************
*Part 1. data import and data manipulation;

filename in 'C:\Users\naz003\work\personal staff\New folder\survival analysis\project.txt';
DATA project;
INFILE in firstobs=2;
INPUT ID VISIT VISITdate yymmdd10. Arm $  CD4 Age parity married HIVpos Contacts Condom;
RUN;

proc sort data=project;by id visit; run;


*Transpose Data From Long Form to the Wide Form**;
data mpl;
	set project;
	by id;
	if first.id;
	keep id Arm cd4 CD4 age parity married;
run;

%macro trans(var,no);
	proc transpose data=project prefix=&var. out=pwide&no.;
		by id;
		var &var.;
	run;
%mend trans;

%trans(HIVpos,1);
%trans(Contacts,2);
%trans(Condom,3);
%trans(VISITdate,4);



data pwide;
	merge mpl pwide1(drop=_NAME_) pwide2(drop=_NAME_) pwide3(drop=_NAME_) pwide4(drop=_NAME_);
	by id;
run;

data pwide;
	set pwide;
		hivpos=max(of hivpos1-hivpos12);
		lastvisit=max(of visitdate1-visitdate12);
		selavisit=largest(2,of visitdate1-visitdate12);

		time=lastvisit-visitdate1+1;
		contactsmean=mean(of contacts1-contacts12);
		condommean=mean(of condom1-condom12);
		week=time/7;
		year=time/365.25;
	select (arm);
		when ('Jadelle') 	armn=1;
		when ('IUD')		armn=2;
		when ('DMPA')		armn=3;
	end;

run;


data pwide;	
	set pwide;
	**Interval Censor****;
	if hivpos=1 then do;
		upper=time;
		lower=selavisit-visitdate1+1;
	end;
	**Righ Censor*******;
	if hivpos=0 then do;
		upper=.;
		lower=time;
	end;
	logpar=log(parity); 
	parsqu=parity**2;
	parcub=parity**3;
	logp2=log(parsqu);
	if married=10 then married=.;
	if parity=21 then parity=.;
run;



proc format;
	value armn 1='Jadelle' 2='IUD' 3='DMPA';
	value married	1='Married'	0='Otherwise';
	value hivpos	1='HIV Positive by the end of the study'	0='HIV Negative by the end of the study';
run;

data pwide;
	set pwide;
	format married married.
		   hivpos	hivpos.
		   lastvisit MMDDYY10. 
		   selavisit MMDDYY10.;

	label armn='Study Arm A Subject Has Been Randomized to'
		  married='Marital Status at Baseline'
		  Age='Age at Baseline'
		  Parity='Number of Times at Baseline that woman had gien birth previously'
		  HIVPOS='HIV Status by the end of the study'
		  CD4='CD4 T Cell count at Baseline(cells/mm^3)'
;

run;


data baseline;
	set pwide;
	keep id arm cd4 age married parity time hivpos;
run;
***************************************************************************
***************************************************************************
*survival analysis part;
%let output=C:\Users\naz003\work\personal staff\New folder\survival analysis;

ods rtf file="&output.\Final_Descriptive_&sysdate..rtf" ;
proc sort data=pwide;by arm;run;

proc means data=pwide;
	var CD4 Age Parity;
run;

proc means data=pwide;
	var CD4 Age Parity;
	by arm;
run;
proc freq data=pwide;
	table married*arm/nopercent norow cmh;
run;

proc freq data=pwide;
	table arm*hivpos/nopercent nocol cmh;
run;

ods rtf close;




***(iii)Kaplan Meier estimation************;
ods rtf file="&output.\Final_AIM1_univariate1&sysdate..rtf" ;
proc lifetest data=pwide plots=(s ls lls);
	time time*hivpos(0);
	strata arm;
run;

/*****Diagnostics for AFT model****/
%macro test(distr);
proc lifereg data=pwide ;
class arm;
model time*hivpos(0)= arm / distribution =&distr;
PROBPLOT;
output out=b cres=rCi;
run;
%mend test;
%test(exponential);
%test(weibull);
%test(lognormal);

***(iv) Cox proportional hazard model with time-dependent variables and model selection****;

ods rtf file="&output.\Final_AIM1_Multivariate&sysdate..rtf" ;

proc phreg data=pwide ;
class arm married parity;
model time*hivpos(0)= arm age cd4 parity married  condoms contact	/
   rl ties=efron selection=stepwise slentry=0.25 slstay=0.15;
hazardratios arm;

array visitdate{*} visitdate1-visitdate12;
array contacts{*} contacts1-contacts12;
array condom{*}	condom1-condom12;
do j=1 to 12;
if lastvisit>visitdate[j] and visitdate[j] ne . then do; 
contact=contacts[j];
condoms=condom[j];
end;
end;
run;
ods rtf close;

***Test Linearity ****; 
ods rtf file="&output.\Final_AIM1_Multivariate&sysdate..rtf" ;

proc phreg data=pwide ;
class arm married  ;
model time*hivpos(0)= arm  cd4  married /ties=efron rl  ;
hazardratios arm;
output out=temp resmart=rr;
run;

proc gplot data=temp;
	plot rr*cd4;
	symbol i=sm60s v=j font=special;
	label rr='Martingale Residual';
run;
ods rtf close;

****Test whether stratifiction is a better choice;*****************************************************************;

data test;
	input  arm $ married cd4;
	datalines;

	Jadelle	1 705.3
	Jadelle	0 705.3
	DMPA	1 705.3
	DMPA	0 705.3
	IUD		1 705.3
	IUD		0 705.3
;
run;

ods rtf file="&output.\Final_AIM1_Final_Multivariate&sysdate..rtf" ;

proc phreg data=pwide ;
	strata arm;
	class arm;
	model time*hivpos(0)= cd4 married /ties=efron  ;
	baseline out=B2 covariates=test survival=s logsurv=ls loglogs=lls cumhaz=chaz;
run;


symbol1 c=default  		v=none i=steplj l=1;
symbol2 c=blue			 v=none i=steplj l=3;
symbol3 c=red 			v=none i=steplj l=5;
symbol4 c=green 		v=none i=steplj l=2;
symbol5 c=yellow 		v=none i=steplj l=4;
symbol6 c=purple		v=none i=steplj l=6;

proc format;
	value strata 1='Married, Jadelle'
				2='Unmarrid, Jadelle'
				3='Married, DMPA'
				4='Unmarrid, DMPA'
				5='Married, IUD'
				6='Unmarrid, IUD'

			;
data b3;
set b2;
if arm='Jadelle' and married=1 then strata=1;
if arm='Jadelle' and married=0 then strata=2;
if arm='DMPA' and married=1 then strata=3;
if arm='DMPA' and married=0 then strata=4;
if arm='IUD	' and married=1 then strata=5;
if arm='IUD	' and married=0 then strata=6;
format strata strata.;
run;



proc sort data=b3; by strata s;run;


title "Estimated Survival Function using Strata Variables";
proc gplot data=b3;
	plot s*time=strata;
run;

title "Estimated Log Cumulative Hazard Function using Strata Variable";
proc gplot data=b3;
	plot lls*time=strata;
run;
ods rtf close;
