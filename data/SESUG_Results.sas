data resp;
set Data.resp;
SexNum = (Sex = "M");
UniqID = Center*100 + ID;
Trt = (Treatment = "A");
run;

%Adj_WinRatio(	RESP, 
				OUT,
				UniqID,
				Visit1 Visit2 Visit3 Visit4,
				Trt);

%Adj_WinRatio(	RESP, 
				OUT,
				UniqID,
				Visit1 Visit2 Visit3 Visit4,
				Trt,
				BASELINE = Baseline,
				COVARS = Age SexNum,
				STRATA = Center);

%Adj_WinOdds(	RESP, 
				OUT,
				UniqID,
				Visit1 Visit2 Visit3 Visit4,
				Trt,
				BASELINE = Baseline,
				COVARS = Age SexNum,
				STRATA = Center);
				
data skin;
set Data.skin;
Stage4 = (STAGE = 4);
Stage5 = (STAGE = 5);
ID = _N_;
if INV = 5 then center = 1;
else if INV = 6 then center = 2;
else if INV = 8 then center = 3;
else if INV = 9 then center = 4;
else if INV = 10 then center = 5;
else center = 6;
center2 = center;
if center = 4 then center2 = 3;
run;

%Adj_WinRatio(	SKIN, 
				OUT,
				ID,
				R1 R2 R3,
				Trt,
				COVARS = Stage,
				STRATA = Center2);

%Adj_WinOdds(	SKIN, 
				OUT,
				ID,
				R1 R2 R3,
				Trt);

data out_sanon;
	set out;
	SE_WP = SE_logWO * WP * (1-WP);
	Chi_Square_WP = ((WP-0.5) / SE_WP)**2;
	p_WP = 1-probchi(Chi_Square_WP, 1);
run;
