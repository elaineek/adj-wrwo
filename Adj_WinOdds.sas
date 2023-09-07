*************************************************************************************************************
*************************************************************************************************************
*************************************************************************************************************

Adj_WinOdds

Ann Marie Weideman, M.S.
Elaine Kearney Kowalewski, M.S.
Gary G. Koch, Ph.D.

2023-09-07	  	Version 1

Arguments:

dsnin 			SAS dataset containing the analysis data. Must be in wide 
				format such that a participant's repeated responses are in 
				a single row and each response is a separate column.
 
dsnout			Name for output dataset.
 
pid             Variable corresponding to unique participant ID.
 
outcomes        List of the variables (each separated by a space) corresponding 
				to outcomes measured at each visit. The outcomes must have at 
				least an ordinal measurement scale with larger values being 
				better than smaller values. Thus, the outcome can be ordered 
				categories or continuous measurements or dichotomies such as 0 
				or 1 or "no" or "yes".
 
arm             Variable for treatment arm.  Required to be a positive
                integer such that the test treatment arm is ALWAYS higher
                in value than the control arm.
 
baseline        Variable corresponding to outcome measurement at baseline. 
				If not specified, no baseline adjustment is employed (which is default).
 
covars          List of the variables corresponding to the covariates
                (measured at baseline) to be used for adjustment.  These
                covariates must be numeric, and can be measured on a binary,
                categorical, ordered categorical, or continuous scale. If
                not specified, no covariate adjustment is employed (which is default).
 
strata          Variable used for stratification.  If not specified, no
                stratification is utilized (which is default).
 
method          SMALL or LARGE used to denote the method employed. The small 
				sample size method is recommended unless within-stratum sample 
				size is reasonably large (e.g., >= 50), number of visits is small 
				(e.g., <=6), and number of covariates is small (e.g., <=4). 
				If not specified, the default is SMALL.
 
debug           0 does not print analysis details to the log and 1 prints analysis 
				details to the log. If not specified, the default is 0.

*************************************************************************************************************
*************************************************************************************************************
************************************************************************************************************;

%MACRO  Adj_WinOdds(DSNIN,
				DSNOUT,
				PID,
				OUTCOMES,
				ARM,
				BASELINE = NONE,
				COVARS = NONE,
				STRATA = NONE,
				METHOD = SMALL,
				DEBUG = 0); 

	%global s r s_two_r_1 r_1 num_visits h;

	%if %length(&DSNIN) = 0 or %length(&DSNOUT) = 0 %then %do;
		%put ERROR: DSNIN and DSNOUT must be specified;
		%goto stopmac;
	%end;

	%if %length(&OUTCOMES) = 0 %then %do;
		%put ERROR: Must have at least one outcome measurement;
		%goto stopmac;
	%end;

	%if %length(&ARM) = 0 %then %do;
		%put ERROR: Must have treatment arm specified;
		%goto stopmac;
	%end;

	%if ^%index(SMALL LARGE, %upcase(&METHOD)) %then %do;
		%put ERROR: METHOD must be one of SMALL LARGE;
		%goto stopmac;
	%end;

	data _null_;
		 %let dsid=%sysfunc(open(&DSNIN));
		 %let checkPID=%sysfunc(varnum(&dsid,&PID));
		 %let checkARM=%sysfunc(varnum(&dsid,&ARM));
		 %if (&checkPID=0) %then 
		 	%let missing_var = &PID;
		 %else %if (&checkARM=0) %then
		 	%let missing_var = &ARM;
		 %else
		 	%let missing_var = "_NONE";
	run;

	%if &missing_var ne "_NONE" %then %do;
		%put ERROR: &missing_var not in &DSNIN;
		%goto stopmac;
	%end;

	proc freq data = &DSNIN.;
		tables &PID. / noprint out = _idlist;
	run;

	proc sql noprint;
		select max(count) into :maxCount
		from _idlist;
	quit;

	proc datasets nolist;
		delete _idlist;
	quit;

	%if &maxCount >= 2 %then %do;
		%put ERROR: PID must be unique for each observation;
		%goto stopmac;
	%end;

	data _null_;
		call symput("num_visits", left(countw(compbl("&OUTCOMES"), " ")));
		%if &COVARS = NONE %then %do;
			%let num_covars = 0;
		%end;
		%else %do;
			call symput("num_covars", left(countw(compbl("&COVARS"), " ")));
		%end;
	run;


	data _null_;
		%let missing_var2 = "_BLANK";
		%let dsid=%sysfunc(open(&DSNIN));
		%do vis = 1 %to &num_visits;
			%let var = %SCAN("&OUTCOMES.",&vis.," ");
			%if (%sysfunc(varnum(&dsid,&var))=0) %then %do;
		 		%let missing_var2 = &var;
				%goto stopdata;
			%end;
		%end;
		%do vis = 1 %to &num_covars;
			%let var = %SCAN("&COVARS.",&vis.," ");
			%if (%sysfunc(varnum(&dsid,&var))=0) %then %do;
		 		%let missing_var2 = &var;
				%goto stopdata;
			%end;
		%end;
		%if &BASELINE. ne NONE %then %do;
			%if (%sysfunc(varnum(&dsid,&BASELINE))=0) %then %do;
		 		%let missing_var2 = &BASELINE;
				%goto stopdata;
			%end;
		%end;
		%if &STRATA. ne NONE %then %do;
			%if (%sysfunc(varnum(&dsid,&STRATA))=0) %then %do;
		 		%let missing_var2 = &STRATA;
				%goto stopdata;
			%end;
		%end;
		%stopdata:;
	run;

	%if &missing_var2 ne "_BLANK" %then %do;
		%put ERROR: &missing_var2 not in &DSNIN;
		%goto stopmac;
	%end;	

	data _local_dsnin;
		set &DSNIN.;
		%if &STRATA. = NONE %then %do;
			NONE = 1;
		%end;
	run;

	proc freq nlevels data = _local_dsnin;
    	ods exclude all;
      	ods output NLevels = _ntrt;
      	tables &ARM.;
   	run;   
   	ods select all;

   	data _null_;
      	set _ntrt;
      	call symput("num_trt", left(nlevels)); 
   	run; 

   	proc datasets nolist;
		delete _ntrt;
	quit;

   	%if &num_trt ne 2 %then %do;
	   	proc datasets nolist;
			delete _local_dsnin;
		quit;
      	%put ERROR: ARM must have TWO levels;
      	%goto stopmac;
   	%end;   	

   	proc sort data = _local_dsnin out = _local_dsnin_sort;
   		by &ARM.;
   	run;

	proc freq noprint data = _local_dsnin_sort;
    	tables &ARM. / out = _ntrts;
	run;

	data _null_;
		set _ntrts;
		if _n_ = 1 then do;
			call symput("ctrl_arm", trim(left(&ARM.)));
			call symput("n_ctrl_arm", trim(left(COUNT)));
		end;
		else if _n_ = 2 then do;
			call symput("trt_arm", trim(left(&ARM.)));
			call symput("n_trt_arm", trim(left(COUNT)));
		end;
	run;

   	proc datasets nolist;
		delete _ntrts;
	quit;

   	proc sql noprint;
      	select distinct &STRATA. 
      	into :stratalist separated by ','
      	from _local_dsnin;
   	quit; 

	proc freq nlevels data = _local_dsnin;
    	ods exclude all;
      	ods output NLevels = _nstrt;
      	tables &STRATA.;
   	run;   
   	ods select all;

   	data _null_;
      	set _nstrt;
      	call symput("num_strata", left(nlevels)); 
   	run; 

   	proc datasets nolist;
		delete _nstrt;
	quit;

   	%if &num_strata. > 0 & &METHOD. = LARGE %then %do;

   		proc freq data = _local_dsnin noprint;
   			tables &STRATA. / out = _nnstrata;
   		run;

		proc means data = _nnstrata noprint;
			var COUNT;
			output out = _min_strata Min = min;
		run;

	   	proc datasets nolist;
			delete _nnstrata;
		quit;

		data _null_;
			set _min_strata;
			call symput("min_ss_strata", left(min));
		run;

	   	proc datasets nolist;
			delete _min_strata;
		quit;

		%if &min_ss_strata. < 50 %then %do;
      		%put WARNING: Minimum within-stratum sample size is 50. Consider using METHOD = SMALL instead.;
   		%end;  

   	%end;

	%if &num_strata. > 0 & &BASELINE. ^= NONE %then %do;

		proc means data = _local_dsnin noprint;
			class &STRATA. &ARM.;
			var &BASELINE.;
			output out = _mean_strata_base stddev = sd;
		run;

		proc means data = _mean_strata_base noprint;
			where _TYPE_ = 3;
			var sd;
			output out = _min_strata_base_sd Min = min;
		run;

	   	proc datasets nolist;
			delete _mean_strata_base;
		quit;

		data _null_;
			set _min_strata_base_sd;
			call symput("min_sd_base_strata", left(min));
		run;

	   	proc datasets nolist;
			delete _min_strata_base_sd;
		quit;

		%if &min_sd_base_strata. = 0 %then %do;
		   	proc datasets nolist;
				delete _local_dsnin _local_dsnin_sort;
			quit;
      		%put ERROR: A treatment group in a strata has the same value for all patients for &BASELINE.;
			%goto stopmac;
   		%end;  

	%end;

	%if &num_strata. > 0 %then %do;

		%do vis = 1 %to &num_visits.;

			proc means data = _local_dsnin noprint;
				class &STRATA. &ARM.;
				var %SCAN("&OUTCOMES.",&vis.," ");
				output out = _mean_strata_out stddev = sd;
			run;

			proc means data = _mean_strata_out noprint;
				where _TYPE_ = 3;
				var sd;
				output out = _min_strata_out_sd Min = min;
			run;

		   	proc datasets nolist;
				delete _mean_strata_out;
			quit;

			data _null_;
				set _min_strata_out_sd;
				call symput("min_sd_out_strata", left(min));
			run;

		   	proc datasets nolist;
				delete _min_strata_out_sd;
			quit;

			%if &min_sd_out_strata. = 0 %then %do;
			   	proc datasets nolist;
					delete _local_dsnin _local_dsnin_sort;
				quit;
	      		%put ERROR: A treatment group in a strata has the same value for all patients for %SCAN("&OUTCOMES.",&vis.," ").;
				%goto stopmac;
	   		%end;  

		%end;

	%end;

	proc sort data = _local_dsnin out = _local_dsnin_sort;
		by &ARM. &STRATA.;
	run;

	proc freq data = _local_dsnin_sort noprint;
		tables &ARM. * &STRATA. / out = _ntrtsstrata;
	run;

	data _null_;
		set _ntrtsstrata;
		if &ARM. = "&ctrl_arm." then do;
			%do i = 1 %to &num_strata.;
				if &STRATA. = %SCAN("&stratalist.",&i.,",") then 
					call symput("ctrl_arm_strata&i.", trim(left(COUNT)));
			%end;
		end;
		else if &ARM. = "&trt_arm." then do;
			%do i = 1 %to &num_strata.;
				if &STRATA. = %SCAN("&stratalist.",&i.,",") then 
					call symput("trt_arm_strata&i.", trim(left(COUNT)));
			%end;
		end;
	run;

	proc datasets nolist;
		delete _ntrtsstrata _local_dsnin_sort;
	quit;

	%if &debug. %then %do;
		%PUT num visits = &num_visits.;
		%PUT num covariates = &num_covars.;
		%PUT control arm = &ctrl_arm.;
		%PUT test arm = &trt_arm.;
		%PUT treatment counts = &n_ctrl_arm. &n_trt_arm.;
		%PUT strata: &stratalist.;
		%PUT num strata = &num_strata.;
		%do h = 1 %to &num_strata.;
			%put num control arm strata &h. = &&ctrl_arm_strata&h;
			%put num test arm strata &h. = &&trt_arm_strata&h;
		%end;
	%end;

	%if &BASELINE = NONE %then %do;
		%let r_1 = &num_visits.;
	%end;
	%else %do;
		%let r_1 = %eval(&num_visits.+1);
	%end;
	%let r = &num_visits.;
	%let s = &num_covars.;
	%let s_two_r_1 = %eval(&s. + 2*&r_1.);

	PROC IML;

		/******************************* COMPUTEU_ODDS FUNCTION *********************************/
		start computeU_odds;
			U_iip = J(ss1*ss2, &s_two_r_1., 0);
			U_i   = J(ss2, &s_two_r_1., 0);
			U_ip  = J(ss1, &s_two_r_1., 0);
			U_h   = J(1, &s_two_r_1., 0);

			counter = 1;
			do i = 1 to ss2;
				do ip = 1 to ss1;

					loc = 1;
					%if &s. > 0 %then %do;
						do k = 1 to &s.;
							U_iip[counter, loc] = trtDataCovar[i, k] - ctrlDataCovar[ip, k];
							loc = loc + 1;
						end;
					%end;

					%if &r_1. > &r. %then %do;
						if (trtData0[i] = . | ctrlData0[ip] = . | (trtData0[i] = ctrlData0[ip])) then 
							U_iip[counter, loc] = 0.5;
						else if (trtData0[i] > ctrlData0[ip]) then
							U_iip[counter, loc] = 1;
						else
							U_iip[counter, loc] = 0;
						loc = loc + 1;
					%end;

					%do j = 1 %to &num_visits.;

						if (trtData[i, &j.] = . | ctrlData[ip, &j.] = . | (trtData[i, &j.] = ctrlData[ip, &j.])) then
							U_iip[counter, loc] = 0.5;
						else if (trtData[i, &j.] > ctrlData[ip, &j.]) then
							U_iip[counter, loc] = 1;
						else
							U_iip[counter, loc] = 0;
						loc = loc + 1;
		
					%end;

					%if &r_1. > &r. %then %do;
						if (trtData0[i] = . | ctrlData0[ip] = . | (trtData0[i] = ctrlData0[ip])) then 
							U_iip[counter, loc] = 0.5;
						else if (trtData0[i] < ctrlData0[ip]) then
							U_iip[counter, loc] = 1;
						else 
							U_iip[counter, loc] = 0;
						loc = loc + 1;
					%end;

					%do j = 1 %to &num_visits.;

						if (trtData[i, &j.] = . | ctrlData[ip, &j.] = . | (trtData[i, &j.] = ctrlData[ip, &j.])) then 
							U_iip[counter, loc] = 0.5;
						else if (trtData[i, &j.] < ctrlData[ip, &j.]) then
							U_iip[counter, loc] = 1;
						else
							U_iip[counter, loc] = 0;
						loc = loc + 1;
		
					%end;

					U_i[i, ] = U_i[i, ] + U_iip[counter, ];
					U_ip[ip, ] = U_ip[ip, ] + U_iip[counter, ];

					U_h = U_h + U_iip[counter, ];

					counter = counter + 1;
				end;
			end;
			U_i_dot   = U_i / ss1;
			U_dot_ip  = U_ip / ss2;

			U_h       = U_h / (ss1*ss2);
			if &debug. then
				print U_h;

			U_i_dot_term = t(U_i_dot - U_h) * (U_i_dot - U_h);
			if &debug. then
				print U_i_dot_term;
			U_dot_ip_term = t(U_dot_ip - U_h) * (U_dot_ip - U_h);
			V_h = (1/(ss2*(ss2-1))) * U_i_dot_term
			  + (1/(ss1*(ss1-1))) * U_dot_ip_term;
			if &debug. then
				print V_h;
		finish computeU_odds;
		/*****************************************************************************************/

		U 	= J(1, &s_two_r_1., 0);
		V 	= J(&s_two_r_1., &s_two_r_1., 0);
		b 	= J(&r., 1, 0);
		V_b = J(&r., &r., 0);

		w_h = J(&num_strata., 1, 0);
		%do h = 1 %to &num_strata.;
			w_h[&h] = (&&ctrl_arm_strata&h*&&trt_arm_strata&h)/(&&ctrl_arm_strata&h+&&trt_arm_strata&h+1);
		%end;
		sum_w_h = w_h[+];
		w_h = w_h / sum_w_h;
		if &debug. then
			print w_h;

		%do h = 1 %to &num_strata.;

			* read covariates *;
			%if &s. > 0 %then %do;
				use _local_dsnin;
					read all var {&COVARS.}
					into ctrlDataCovar 
					where (&ARM. = &ctrl_arm. & &STRATA. = %SCAN("&stratalist.",&h.,","));
				close _local_dsnin;

				use _local_dsnin;
					read all var {&COVARS.}
					into trtDataCovar
					where (&ARM. = &trt_arm. & &STRATA. = %SCAN("&stratalist.",&h.,","));
				close _local_dsnin;
			%end;
			* read baseline *;
			%if &r_1. > &r. %then %do;
				use _local_dsnin;
					read all var {&BASELINE.}
					into ctrlData0
					where (&ARM. = &ctrl_arm. & &STRATA. = %SCAN("&stratalist.",&h.,","));
				close _local_dsnin;

				use _local_dsnin;
					read all var {&BASELINE.}
					into trtData0 
					where (&ARM. = &trt_arm. & &STRATA. = %SCAN("&stratalist.",&h.,","));
				close _local_dsnin;
			%end;

			use _local_dsnin;
				read all var {&OUTCOMES.}
				into ctrlData
				where (&ARM. = &ctrl_arm. & &STRATA. = %SCAN("&stratalist.",&h.,","));
			close _local_dsnin;

			use _local_dsnin;
				read all var {&OUTCOMES.}
				into trtData 
				where (&ARM. = &trt_arm. & &STRATA. = %SCAN("&stratalist.",&h.,","));
			close _local_dsnin;

			ss1 = &&ctrl_arm_strata&h;
			ss2 = &&trt_arm_strata&h;

			RUN computeU_odds;

			%if &METHOD. = SMALL %then %do;
				U = U + w_h[&h.] * U_h;
				V = V + w_h[&h.]*w_h[&h.] * V_h;
			%end;
			%else %do;
				U_h = t(U_h);
				Aright = -1 * I(&r_1.);
				A = I(&r_1.) || Aright;
				U_12 = U_h[(&s.+1):nrow(U_h)];
				if &debug. then
					print U_12;
				f = A * log (U_12);
				%if &r_1. > &r. %then %do;
					f_0 = f[1];
					f_star = f[2:nrow(f)];
				%end;
				%else %do;
					f_star = f[1:nrow(f)];
				%end;
				D = diag(U_12);
				ADinv = A * inv(D);
				%if &s. > 0 %then %do;
					U_x = U_h[1:&s.];
					capF = U_x // f;
					L1top = I(&s.) || J(&s., 2*&r_1., 0);
					L1bot = J(&r_1., &s., 0) || ADinv;
					L1 = L1top // L1bot;
				%end;
				%else %do;
					capF = f;
					L1 = ADinv;
				%end;
				%if &r_1. > &r. %then %do;
					L = J(&s.+1, &r., 0) // I(&r.);
				%end;
				%else %if &s. > 0 %then %do;
					L = J(&s., &r., 0) // I(&r.);
				%end;
				%else %do; * no covariates and no baseline *;
					L = I(&r.);
				%end;
				V_F = L1 * V_h * t(L1);
				b_h = inv(t(L) * inv(V_F) * L) * t(L) * inv(V_F) * capF;
				V_b_h = inv(t(L) * inv(V_F) * L);

				b = b + w_h[&h.] * b_h;
				V_b = V_b + w_h[&h.]*w_h[&h.] * V_b_h;
			%end;
		%end;

		%if &METHOD. = SMALL %then %do;
			U = t(U);
			if &debug. then
				print U;
			if &debug. then
				print V;

			Aright = -1 * I(&r_1.);
			A = I(&r_1.) || Aright;
			U_12 = U[(&s.+1):nrow(U)];
			f = A * log (U_12);
			if &debug. then
				print f;
			%if &r_1. > &r. %then %do;
				f_0 = f[1];
				f_star = f[2:nrow(f)];
			%end;
			%else %do;
				f_star = f[1:nrow(f)];
			%end;
			D = diag(U_12);
			ADinv = A * inv(D);
			%if &s. > 0 %then %do;
				U_x = U[1:&s.];
				capF = U_x // f;
				L1top = I(&s.) || J(&s., 2*&r_1., 0);
				L1bot = J(&r_1., &s., 0) || ADinv;
				L1 = L1top // L1bot;
			%end;
			%else %do;
				capF = f;
				L1 = ADinv;
			%end;
			%if &r_1. > &r. %then %do;
				L = J(&s.+1, &r., 0) // I(&r.);
			%end;
			%else %if &s. > 0 %then %do;
				L = J(&s., &r., 0) // I(&r.);
			%end;
			%else %do; * no covariates and no baseline *;
				L = I(&r.);
			%end;
			V_F = L1 * V * t(L1);

			if &debug. then
				print V_F;

			b = inv(t(L) * inv(V_F) * L) * t(L) * inv(V_F) * capF;
			V_b = inv(t(L) * inv(V_F) * L);
		%end;

		if &debug. then
			print b;

		if &debug. then
			print V_b;
		
		var_b = vecdiag(V_b);
		se_b = sqrt(var_b);

		chi_b = (b / sqrt(var_b)) ## 2;
		p_b = 1-probchi(chi_b, 1);

		logWOj =  b||se_b||chi_b||p_b;
		create _outds from logWOj[colname={"logWO" "SE logWO" "Chi Square" "p-value"}];
		append from logWOj;

		create &DSNOUT._b from b[colname={'b'}];
		append from b;

		create &DSNOUT._Vb from V_b[colname={&OUTCOMES.}];
		append from V_b;

	QUIT;

	data _Outcomes;
		%do i=1 %to &num_visits.;
		Visit = "%scan(&OUTCOMES., &i)";
		output;
		%end;
	run;
	
	data &DSNOUT;
		merge _Outcomes _outds;
		WO = exp(logWO);
		WO_CI = cats("(", put(exp(logWO - 1.96*SE_logWO), 6.2), ", ", put(exp(logWO + 1.96*SE_logWO), 6.2), ")");
		WP = WO / (1 + WO);
		format logWO 6.3
			   SE_logWO 6.3
			   Chi_Square 6.2
			   p_value pvalue5.3
			   WO 6.2
			   WP 6.3;
		label logWO = "log(WO)"
			  SE_logWO = "SE log(WO)"
			  Chi_Square = "Chi-square"
			  p_value = "P-value"
			  WO_CI = "WO 95% CI";
	run;

	proc print data = &DSNOUT label noobs;
	run;

	proc datasets nolist;
		delete _Outcomes _outds _local_dsnin;
	quit;
	%stopmac:;
%MEND Adj_WinOdds;
