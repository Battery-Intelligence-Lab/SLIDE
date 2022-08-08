// void Cell_SPM::setStates(State_SPM &&states)
// {
// 	/*
// 	 * sets all battery state variables
// 	 *
// 	 * IN
// 	 * states	array with the battery states
// 	 *
// 	 */

// 	if constexpr (settings::printBool::printCellFunctions)
// 		std::cout << "Cell_SPM::setStates(int, double[]) starting.\n";

// 	validState(); // throws if the state is illegal

// 	st = std::move(states);

// 	// the stress values stored in the class variables for stress are no longer valid because the state has changed
// 	sparam.s_dai_update = false;
// 	sparam.s_lares_update = false;

// 	if constexpr (settings::printBool::printCellFunctions)
// 		std::cout << "Cell_SPM::setStates(int, double[]) terminating.\n";
// }

// void Cell_SPM::setStates(const State_SPM &si, double I)
// {
// 	/*
// 	 * Set all battery state variables to the ones in the State-object provided.
// 	 * Set the cell current to the value provided (without ramping the current since we are changing the state any way)
// 	 *
// 	 * IN
// 	 * si 	new state of the cell
// 	 * I 	new current of the cell [A], positive for discharging, negative for charging
// 	 */

// 	if constexpr (settings::printBool::printCellFunctions)
// 		std::cout << "Cell_SPM::setStates(State, double) starting.\n";

// 	st = si;
// 	Icell = I;

// 	// the stress values stored in the class variables for stress are no longer valid because the state has changed
// 	sparam.s_dai_update = false;
// 	sparam.s_lares_update = false;

// 	if constexpr (settings::printBool::printCellFunctions)
// 		std::cout << "Cell_SPM::setStates(State, double) terminating.\n";
// }

void Cell_SPM::peekVoltage(double I)
{
    /*
     * Function to set the cell current to the specified value then peek the voltage
     * The current is slowly ramped from the actual cell current to the specified value.
     *
     * IN
     * I 		value to which the current should be set [A]
     * 			> 0 for discharge
     * 			< 0 for charge
     *
     * THROWS
     * 105 		check == true and the specified current could not be set without violating the cell's limits.
     * 			the original battery state and current are restored
     */
    // if constexpr (settings::printBool::printCellFunctions)
    // 	std::cout << "Cell_SPM::setCurrent starting with current " << I << '\n';

    // // check if the specified value is different from the actual cell current
    // if (std::abs(Icell - I) < 1e-10)
    // 	return; // the values are the same -> we don't need to do anything

    // // Store the old current and state to restore it if needed
    // const auto sold = getStates();
    // const auto Iold = I();

    // // settings
    // bool blockDegradation = true; // don't account for degradation while the current is ramping
    // bool reached = false;		  // boolean indicating if we have reached the set current
    // int sgn;					  // sgn whether we should increase or decrease the current
    // if (I > Icell)
    // 	sgn = 1;
    // else
    // 	sgn = -1;

    // // loop to ramp the current
    // while (!reached)
    // {
    // 	// increase the current
    // 	Icell += sgn * dIcell;

    // 	// check if you have reached the specified current
    // 	if (sgn * Icell > sgn * I)
    // 	{ // increase I: Icell > I, decrease I: Icell < I
    // 		Icell = I;
    // 		reached = true;
    // 	}

    // 	// take one small time step
    // 	ETI(false, dt_I, blockDegradation);
    // }

    // // Check the cell's conditions are still valid if we want to check the final state
    // double v, ocvp, ocvn, etap, etan, rdrop, tem;
    // bool valid;
    // if (check)
    // {
    // 	// check the state
    // 	try
    // 	{
    // 		validState(); // throws an error if the state is illegal
    // 	}
    // 	catch (int e)
    // 	{
    // 		if (print)
    // 			std::cout << "Cell_SPM::setCurrent illegal state after setting the current to " << Icell << ", error: " << e << ". Throwing an error.\n";
    // 		setStates(sold, Iold); // restore the original battery state and current
    // 		throw e;
    // 	}

    // 	// check the voltage
    // 	try
    // 	{
    // 		const int status = getVoltage(print, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
    // 		valid = (status == 0); // throws an error if the surface concentration is out of bounds
    // 	}
    // 	catch (int e)
    // 	{
    // 		valid = false; // in that case, the voltage is illegal
    // 	}
    // 	if (!valid)
    // 	{
    // 		if (print)
    // 			std::cerr << "Cell_SPM::setCurrent Illegal voltage after trying to set the current to " << I() << ", the voltage is: " << v << "V. Throwing an error.\n";
    // 		setStates(sold, Iold); // restore the original battery state and current
    // 		throw 105;
    // 	}
    // }

    // // the stress values stored in the class variables for stress are no longer valid because the cell current has changed
    // sparam.s_dai_update = false;
    // sparam.s_lares_update = false;

    // if constexpr (settings::printBool::printCellFunctions)
    // {
    // 	if (check)
    // 		std::cout << "Cell_SPM::setCurrent terminating with current " << I << " and voltage " << v << '\n';
    // 	else
    // 		std::cout << "Cell_SPM::setCurrent terminating with current " << I << " without checking the voltage.\n";
    // }
}

// void Cell_SPM::getStates(State_SPM &si, double *I)
// {
// 	/*
// 	 * Copies all the states which describe the status of the cell to the State-object
// 	 *
// 	 * OUT
// 	 * si 	a reference to a state-object in which the states will be written
// 	 * I 	the actual cell current [A]
// 	 * 			positive for discharging
// 	 * 			negative for charging
// 	 */

// 	// #HOTFUNC  -> called 450 million times.
// 	if constexpr (settings::printBool::printCellFunctions)
// 		std::cout << "Cell_SPM::getStates(State, double*) starting\n";

// 	si = st;
// 	*I = Icell;

// 	if constexpr (settings::printBool::printCellFunctions)
// 		std::cout << "Cell_SPM::getStates(State, double*) terminating.\n";
// }

// // state space model
// States_SPM Cell_SPM::dState_SLIDEv2(bool print, bool blockDegAndTherm, int electr)
// {
// 	/*
// 	 * function calculating the time derivatives of the battery states.
// 	 *
// 	 * IN
// 	 * print 	boolean indicating if we want to print error messages or not
// 	 * 				if true, error messages are printed
// 	 * 				if false no error messages are printed (but the error will still be thrown)
// 	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
// 	 * blockDegAndTherm 	if true, battery degradation is ignored (i.e. the time derivatives of those states are 0)
// 	 * electr 	integer describing which electrodes should be considered in the case you want to do half-cell cycling.
// 	 * 			half-cell cycling is only allowed if blockDegAndTherm is true (i.e. you ignore degradation)
// 	 * 			this is because some of the degradation mechanisms could give NaN or Inf due to a divide by 0
// 	 * 				1 					only the positive electrode is considered [half-cell cycling with only the positive electrode]
// 	 * 				2 					only the negative electrode is considered [half-cell cycling with only the negative electrode]
// 	 * 				any other value		both electrodes are considered [normal operation]
// 	 *
// 	 * OUT
// 	 * dstates	change in the states
// 	 * 		dzp			time derivative of the transformed concentration at the positive inner nodes of the positive electrode (dzp/dt)
// 	 * 		dzn			time derivative of the transformed concentration at the positive inner nodes of the negative electrode (dzn/dt)
// 	 * 		dT			time derivative of the battery temperature [K s-1] (dT/dt)
// 	 * 		ddelta 		time derivative of the SEI thickness [m s-1] (ddelta/dt)
// 	 * 		dLLI 		time derivative of the lost lithium inventory [C s-1] (dLLI/dt)
// 	 * 		dthickp 	time derivative of the thickness of the positive electrode [m s-1] (dthickp/dt), <0 (dthickp/dt)
// 	 * 		dthickn		time derivative of the thickness of the negative electrode [m s-1] (dthickn/dt), <0 (dthickn/dt)
// 	 * 		dep			time derivative of the volume fraction in the positive electrode [s-1] (dep/dt)
// 	 * 		den			time derivative of the volume fraction in the negative electrode [s-1] (den/dt)
// 	 * 		dap			time derivative of the effective surface area of the positive electrode [m2 m-3 s-1] (dap/dt)
// 	 * 		dan			time derivative of the effective surface area of the negative electrode [m2 m-3 s-1] (dan/dt)
// 	 * 		dCS 		time derivative of the crack surface [m2 s-1], dCS/st > 0 (dCS/dt)
// 	 * 		dDp 		time derivative of the diffusion constant at reference temperature of the positive electrode [m s-1 s-1] (dDp/dt)
// 	 * 		dDn			time derivative of the diffusion constant at reference temperature of the negative electrode [m s-1 s-1] (dDn/dt)
// 	 * 		dR			time derivative of the electrode resistance [Ohm m2 s-1] (dR/dt)
// 	 * 		ddelta_pl 	time derivative of the thickness of the plated lithium layer [m s-1] (ddelta_pl/dt)
// 	 *
// 	 * THROWS
// 	 * 101 	the surface li-concentration is out of bounds
// 	 * 109	illegal (combination of) input parameters
// 	 */

// 	using namespace PhyConst;
// 	using settings::nch;

// 	if constexpr (settings::printBool::printCellFunctions)
// 		std::cout << "Cell_SPM::dState starting\n";

// 	// Update the stress values stored in the attributes with the stress of the previous time step
// 	sparam.s_dai_p_prev = sparam.s_dai_p;	  // Dai's stress in the positive particle in the previous time step
// 	sparam.s_dai_n_prev = sparam.s_dai_n;	  // Dai's stress in the negative particle in the previous time step
// 	sparam.s_lares_n_prev = sparam.s_lares_n; // Laresgoiti's stress in the negative particle in the previous time step

// 	// #TODO: updateDaiStress and updateLaresgoitiStress may not need to be calculated if blockDegAndTherm is true.
// 	// Calculate the stress values stored in the attributes for the stress in this time step
// 	if (sparam.s_dai) // only a few degradation models actually need the stress according to Dai, so only calculate it if needed
// 		updateDaiStress();
// 	if (sparam.s_lares) // only a few degradation models need the stress according to Laresgoiti
// 		updateLaresgoitiStress(print);

// 	if ((electr == 1 || electr == 2) && !blockDegAndTherm)
// 	{
// 		std::cerr << "ERROR in Cell_SPM::dState. you are cycling with only one electrode " << electr
// 				  << " but you are also accounting for degradation."
// 					 " That is not allowed, half-cell cycling is only allowed if no degradation is considered. Either change the value of 'electr' to 3"
// 					 " such that you cycle with the full cell, or set 'blockDegAndTherm' to true to ignore degradation.\n";
// 		throw 109;
// 	}

// 	// expanded version of Cell_SPM::getCSurf:

// 	auto [Dpt, Dnt] = calcDiffConstant();
// 	auto [i_app, jp, jn] = calcMolarFlux(); // current density, molar flux on the pos/neg particle
// 	auto [cps, cns] = calcSurfaceConcentration(jp, jn, Dpt, Dnt);

// 	// check if the surface concentration is within the allowed range
// 	// 	0 < cp < Cmaxpos
// 	// 	0 < cn < Cmaxneg
// 	if (cps <= 0 || cns <= 0 || cps >= Cmaxpos || cns >= Cmaxneg) // Do not delete if you cannot ensure zp/zn between 0-1
// 	{
// 		if (print)
// 		{
// 			std::cerr << "ERROR in Cell_SPM::dState: concentration out of bounds. the positive lithium fraction is " << cps / Cmaxpos << " and the negative lithium fraction is " << cns / Cmaxneg;
// 			std::cerr << "they should both be between 0 and 1.\n";
// 		}
// 		throw 101;
// 	}

// 	const double zp_surf = (cps / Cmaxpos); // lithium fraction (0 to 1)
// 	const double zn_surf = (cns / Cmaxneg);
// 	auto [etap, etan] = calcOverPotential(cps, cns, i_app);

// 	if (electr == 2) // only consider negative electrode, ignore the positive electrode
// 	{
// 		jp = 0;
// 		Dpt = 0;
// 		etap = 0;
// 	}
// 	else if (electr == 1) // only consider positive electrode, ignore the negative electrode
// 	{
// 		jn = 0;
// 		Dnt = 0;
// 		etan = 0;
// 	}

// 	// Calculate the entropic coefficient
// 	bool bound = true;													   // in linear interpolation, throw an error if you are outside of the allowed range of the data
// 	const double dOCV = OCV_curves.dOCV_tot.interp(zp_surf, print, bound); // entropic coefficient of the entire cell OCV [V K-1]

// 	// temperature
// 	// Calculate the thermal sources/sinks/transfers per unit of volume of the battery
// 	// The battery volume is given by the product of the cell thickness and the electrode surface
// 	const double Qrev = -i_app / L * st.T() * dOCV;			  // reversible heat due to entropy changes [W m-3]
// 	const double Qrea = i_app / L * (etan - etap);			  // reaction heat due to the kinetics [W m-3]
// 	const double Qohm = I() * I() * getR() / (L * elec_surf); // Ohmic heat due to electrode resistance [W m-3]
// 	const double Qc = -Qch * SAV * (st.T() - T_env);		  // cooling with the environment [W m-3]

// 	// Calculate the effect of the main li-reaction on the (transformed) concentration
// 	State_SPM::z_type dzp, dzn;
// 	for (int j = 0; j < nch; j++)
// 	{											// loop for each row of the matrix-vector product A * z
// 		const double ctep = M.Ap[j] * st.zp(j); // A is diagonal, so the array M.A has only the diagonal elements
// 		const double cten = M.An[j] * st.zn(j);
// 		dzp[j] = (Dpt * ctep + M.Bp[j] * jp); // dz/dt = D * A * z + B * j
// 		dzn[j] = (Dnt * cten + M.Bn[j] * jn);
// 	}

// 	// If we ignore degradation in this time step, we have calculated everything we need
// 	if (blockDegAndTherm)
// 	{
// 		slide::states_type dstates{};									   // Initialize as zero.
// 		std::copy(dzp.begin(), dzp.end(), dstates.begin());				   // first nch dstates are d_zp,
// 		std::copy(dzn.begin(), dzn.end(), dstates.begin() + nch);		   // first nch dstates are d_zn,
// 		dstates[2 * nch + 0] = 1 / (rho * Cp) * (Qrev + Qrea + Qohm + Qc); // dT		cell temperature

// 		// Others are zero : ddelta	SEI thickness, dLLI	lost lithium, dthickp/dthickn 	electrode thickness,
// 		// dep/den volume fraction of active material, dap/dan effective surface are, a = 3 e/R, dCS surface area of the cracks,
// 		// dDp/dDn diffusion constant, dR electrode resistance, ddelta_pl thickness of the plated lithium

// 		if constexpr (settings::printBool::printCellFunctions)
// 			std::cout << "Cell_SPM::dState terminating without degradation.\n";

// 		return dstates; // stop calculating
// 	}

// 	// calculate the anode potential (needed for various degradation models)
// 	const double dOCVn = OCV_curves.dOCV_neg.interp(zn_surf, print, bound); // entropic coefficient of the anode potential [V K-1]
// 	const double OCV_n = OCV_curves.OCV_neg.interp(zn_surf, print, bound);	// anode potential [V]
// 	const double OCVnt = OCV_n + (st.T() - T_ref) * dOCVn;					// anode potential at the cell's temperature [V]

// 	// SEI growth
// 	double isei;		// current density of the SEI growth side reaction [A m-2]
// 	double den_sei;		// decrease in volume fraction due to SEI growth [s-1]
// 	double dznsei[nch]; // additional diffusion in the anode due to isei

// 	SEI(OCVnt, etan, &isei, &den_sei); // Throws but not wrapped in try-catch since only appears here.

// 	// Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF)
// 	for (int j = 0; j < nch; j++)
// 		dznsei[j] = (M.Bn[j] * isei / (nsei * F));

// 	// crack growth leading to additional exposed surface area
// 	double isei_multiplyer; // relative increase in isei due to additional SEI growth on the extra exposed surface area [-]
// 	double dCS;				// increase in crack surface area [m2 s-1]
// 	double dDn;				// change in negative diffusion constant [m s-1 s-1]
// 	double dznsei_CS[nch];	// additional diffusion in the anode due to extra SEI growth on the crack surface

// 	CS(OCVnt, etan, &isei_multiplyer, &dCS, &dDn); // Throws but not wrapped in try-catch since only appears here.

// 	// crack surface leads to extra SEI growth because the exposed surface area increases.
// 	// (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF + isei_CS/nF)
// 	const double isei_CS = isei * isei_multiplyer; // extra SEI side reaction current density due to the crack surface area [A m-2]
// 	for (int j = 0; j < nch; j++)
// 		dznsei_CS[j] = (M.Bn[j] * isei_CS / (nsei * F));

// 	// loss of active material LAM
// 	double dthickp, dthickn, dap, dan, dep, den; // change in geometric parameters describing the amount of active material

// 	LAM(print, zp_surf, etap, &dthickp, &dthickn, &dap, &dan, &dep, &den); // Throws but not wrapped in try-catch since only appears here.

// 	// lithium plating
// 	double ipl;			// current density of the plating side reaction [A m-2]
// 	double dzn_pl[nch]; // additional diffusion in the anode due to ipl
// 	try
// 	{
// 		LiPlating(OCVnt, etan, &ipl);
// 	}
// 	catch (int e)
// 	{
// 		// std::cout << "Throw test: " << 41 << '\n';
// 		if (print)
// 			std::cout << "Error in Cell_SPM::dState when calculating the lithium plating: " << e << ".\n";
// 		throw e;
// 	}

// 	// Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + ipl/nF)
// 	for (int j = 0; j < nch; j++)
// 		dzn_pl[j] = (M.Bn[j] * ipl / (npl * F));

// 	// time derivatives
// 	slide::states_type dstates{}; // Initialize as zero.
// 	for (int j = 0; j < nch; j++)
// 	{
// 		dstates[j] = dzp[j];												// dzp 		diffusion
// 		dstates[nch + j] = (dzn[j] + dznsei[j] + dznsei_CS[j] + dzn_pl[j]); // dzn		jtot = jn + isei/nF + isei_CS/nF + ipl/nF
// 	}
// 	dstates[2 * nch + 0] = 1 / (rho * Cp) * (Qrev + Qrea + Qohm + Qc);				   // dT 		cell temperature
// 	dstates[2 * nch + 1] = isei / (nsei * F * rhosei);								   // ddelta	thickness of the SEI layer
// 																					   // delta uses only isei (and not isei + isei_CS) since crack growth increases the area, not the thickness
// 	dstates[2 * nch + 2] = (isei + isei_CS + ipl) * elec_surf * st.thickn() * st.an(); // dLLI 	loss of lithium
// 																					   // i_sei = density => * active surface area = * (surf*thick*specific_surf_neg)
// 	dstates[2 * nch + 3] = dthickp;													   // dthickp 	electrode thickness
// 	dstates[2 * nch + 4] = dthickn;													   // dthickn
// 	dstates[2 * nch + 5] = dep;														   // dep		volume fraction of active material
// 	dstates[2 * nch + 6] = den + den_sei;											   // den
// 	dstates[2 * nch + 7] = dap + 3 / Rp * dep;										   // dap		effective surface area, a = 3 e/R -> da/dt = da/dt + 3/R de/dt
// 	dstates[2 * nch + 8] = dan + 3 / Rn * (den + den_sei);							   // dan
// 	dstates[2 * nch + 9] = dCS;														   // dCS 		surface area of the cracks
// 	dstates[2 * nch + 10] = 0;														   // dDp 		diffusion constant
// 	dstates[2 * nch + 11] = dDn;													   // dDn
// 	dstates[2 * nch + 12] = 0;														   // dR 		specific electrode resistance
// 	dstates[2 * nch + 13] = ipl / (npl * F * rhopl);								   // ddelta_pl thickness of the plated lithium

// 	if constexpr (settings::printBool::printCellFunctions)
// 		std::cout << "Cell_SPM::dState terminating with degradation.\n";
// 	return dstates;
// }

void validState(State_SPM &st, State_SPM &s_ini)
{
    /* Moved from State to Cell.
     * Check if this State object has valid parameters.
     * Errors are thrown if they are not valid.
     *
     * Note to the user: most of these limits are purely intended to indicate something strange is going on.
     * I.e. exceeding the limits below normally doesn't lead to errors in the code.
     * They just indicate things which should not be happening (e.g. increasing the amount of active material).
     * If a limit is 'critical' (i.e. exceeding this limit will give errors in the code), this is indicated in comments below
     *
     * The main purpose of the limits on degradation parameters (e.g. too low amount of active material)
     * is to avoid a very long, or even infinitely long, calculation (e.g. when the cell simply can't reach a certain voltage because there is not enough active material)
     * or to avoid problems with the time discretisation
     * (e.g. the lower amount of active material, the larger the concentration difference for constant current and time step;
     * at some point, the concentration difference over one time step will become too large, and lead to numerical instability or other discretisation errors)
     *
     * THROWS
     * 15 	illegal state. The error message will detail which state has an illegal value
     */

    // There are no limits on the transformed concentration, because this is the (twice) transformed concentration.
    // The limits are on the real concentrations, which must be calculated first.
    // This can't be done here because it requires parameters of the Cell itself.
    // see Cell_SPM::getC or Cell_SPM::getCsurf

    // temperature of the cell (in Kelvin), the temperature limits are defined in State.hpp
    if (st.T() < Tmin()) // check the temperature is above 0 degrees
    {
        std::cerr << "Error in State::setT. The temperature " << st.T()
                  << "K is too low. The minimum value is " << settings::Tmin_K << '\n';
        throw 15;
    }
    else if (st.T() > Tmax()) // check the temperature is below 60 degrees
    {
        std::cerr << "Error in State::setT. The temperature " << st.T()
                  << "K is too high. The maximum value is " << settings::Tmax_K << '\n';
        throw 15;
    }

    // thickness of the SEI layer
    const bool del = st.delta() <= 0;
    if (del)
    {
        std::cerr << "Error in State::validState. The SEI thickness delta is " << st.delta()
                  << ", which is too low. Only strictly positive values are allowed.\n";
        throw 15;
    }
    // a value of 0 gives problems in some equations, which have a term 1/delta, which would give nan or inf if delta = 0
    // a negative value might lead to a further decrease in SEI thickness, so it will keep getting more and more negative

    // lost lithium
    bool li = st.LLI() < 0;
    if (li)
    {
        std::cerr << "Error in State::validState. The lost lithium LLI is " << st.LLI()
                  << ", which is too low. Only non-negative values are allowed.\n";
        throw 15;
    }
    // a value of 0 is allowed (it won't lead to errors in the code)

    // thickness of the electrodes
    bool tpmin = st.thickp() <= s_ini.thickp() / 5;
    if (tpmin)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the thickness of the positive electrode is " << st.thickp() << ", which is too low. The minimum is " << s_ini.thickp() / 5 << ", 1/5 of the original thickness"
                  << '\n';
        throw 15;
    }
    // errors will happen if the thickness becomes 0 or negative. otherwise no direct problems are expected
    // on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
    bool tpmax = st.thickp() > s_ini.thickp() * 1.001; // leave a 0.1% margin for numerical errors
    if (tpmax)
    {
        std::cerr << "Error in State::validState. The thickness of the positive electrode is " << st.thickp()
                  << ", which is too high. The maximum is " << s_ini.thickp() << ", the original thickness\n";
        throw 15;
    }
    bool tnmin = st.thickn() <= s_ini.thickn() / 5;
    if (tnmin)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the thickness of the negative electrode is " << st.thickn()
                  << ", which is too low. The minimum is " << s_ini.thickn() / 5 << ", 1/5 of the original thickness\n";
        throw 15;
    }

    // errors will happen if the thickness becomes 0 or negative. otherwise no direct problems are expected
    // on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
    bool tnmax = st.thickn() > s_ini.thickn() * 1.001;
    if (tnmax)
    {
        std::cerr << "Error in State::validState. The thickness of the negative electrode is " << st.thickn()
                  << ", which is too high. The maximum is " << s_ini.thickn() << ", the original thickness.\n";

        throw 15;
    }

    // volume fraction of the active material in the electrodes
    bool epmin = st.ep() <= s_ini.ep() / 5;
    if (epmin)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the volume fraction of the positive electrode is "
                  << st.ep() << ", which is too low. The minimum is " << s_ini.ep() / 5 << ", 1/5 of the original volume fraction.\n";
        throw 15;
    }

    // errors will happen if the volume fraction becomes 0 or negative. otherwise no direct problems are expected
    // on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
    bool epmax = st.ep() > s_ini.ep() * 1.001;
    if (epmax)
    {
        std::cerr << "Error in State::validState. The volume fraction of the positive electrode is " << st.ep()
                  << ", which is too high. The maximum is " << s_ini.ep() << ", the original volume fraction.\n";
        throw 15;
    }

    bool enmin = st.en() <= s_ini.en() / 5;
    if (enmin)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the volume fraction of the negative electrode is " << st.en()
                  << ", which is too low. The minimum is " << s_ini.en() / 5 << ", 1/5 of the original volume fraction.\n";
        throw 15;
    }

    // errors will happen if the volume fraction becomes 0 or negative. otherwise no direct problems are expected
    // on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
    bool enmax = st.en() > s_ini.en() * 1.001;
    if (enmax)
    {
        std::cerr << "Error in State::validState. The volume fraction of the negative electrode is " << st.en()
                  << ", which is too high. The maximum is " << s_ini.en() << ", the original volume fraction.\n";
        throw 15;
    }

    // effective surface area of the porous electrodes
    bool apmin = st.ap() <= s_ini.ap() / 5;
    if (apmin)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the effective surface of the positive electrode is " << st.ap()
                  << ", which is too low. The minimum is " << s_ini.ap() / 5 << ", 1/5 of the original effective surface.\n";
        throw 15;
    }
    // errors will happen if the effective surface becomes 0 or negative. otherwise no direct problems are expected
    // on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
    bool apmax = st.ap() > s_ini.ap() * 1.001;
    if (apmax)
    {
        std::cerr << "Error in State::validState. The effective surface of the positive electrode is " << st.ap()
                  << ", which is too high. The maximum is " << s_ini.ap() << ", the original effective surface.\n";
        throw 15;
    }

    bool anmin = st.an() <= s_ini.an() / 5;
    if (anmin)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the effective surface of the negative electrode is " << st.an()
                  << ", which is too low. The minimum is " << s_ini.an() / 5 << ", 1/5 of the original effective surface.\n";
        throw 15;
    }

    // errors will happen if the effective surface becomes 0 or negative. otherwise no direct problems are expected
    // on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
    bool anmax = st.an() > s_ini.an() * 1.001;
    if (anmax)
    {
        std::cerr << "Error in State::validState. The effective surface of the negative electrode is " << st.an()
                  << ", which is too high. The maximum is " << s_ini.an() << ", the original effective surface.\n";
        throw 15;
    }

    // surface area of the cracks growing at the particle surface
    bool csmin = st.CS() <= 0;
    if (csmin)
    {
        std::cerr << "Error in State::validState. The crack surface area is " << st.CS()
                  << ", which is too low. It must be strictly positive.\n";
        throw 15;
    }

    // don't allow 0 because it might give problems in some of the equations, which might grow CS proportional to the existing CS (so a value of 0 gives 0 increase)
    // a negative value might lead to a further decrease in CS, so it will keep getting more and more negative
    bool csmax = st.CS() > 1e4 * s_ini.CS();
    if (csmax)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the crack surface area is "
                  << st.CS() << ", which is too high. The maximum is " << 1e4 * s_ini.CS()
                  << ", 10,000 times the original crack surface area.\n";
        throw 15;
    }
    // normally, the initial value is 1% of the total real electrode surface area, so 10,000*initial value = 100 * total electrode surface area
    // but in theory no maximum value will give errors in the code

    // diffusion constant at reference temperature for the electrodes
    bool dpmin = st.Dp() <= s_ini.Dp() / 5;
    if (dpmin)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the diffusion constant of the positive electrode is " << st.Dp()
                  << ", which is too low. The minimum is " << s_ini.Dp() / 5 << ", 1/5 of the original diffusion constant.\n";
        throw 15;
    }

    // errors will happen if the diffusion constant becomes 0 or negative. otherwise no direct problems are expected
    // on a longer term, you will get discretisation errors (see above, too bad diffusion -> too large surface concentration swings in one time step -> numerical errors / concentration out of bound, etc
    bool dpmax = st.Dp() > s_ini.Dp() * 1.001;
    if (dpmax)
    {
        std::cerr << "Error in State::validState. The diffusion constant of the positive electrode is "
                  << st.Dp() << ", which is too high. The maximum is "
                  << s_ini.Dp() << ", the original diffusion constant.\n";
        throw 15;
    }

    bool dnmin = st.Dn() <= s_ini.Dn() / 5;
    if (dnmin)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the diffusion constant of the negative electrode is " << st.Dn()
                  << ", which is too low. The minimum is " << s_ini.Dn() / 5 << ", 1/5 of the original diffusion constant.\n";
        throw 15;
    }
    // errors will happen if the diffusion constant becomes 0 or negative. otherwise no direct problems are expected
    // on a longer term, you will get discretisation errors (see above, too bad diffusion -> too large surface concentration swings in one time step -> numerical errors / concentration out of bound, etc
    bool dnmax = st.Dn() > s_ini.Dn() * 1.001;
    if (dnmax)
    {
        std::cerr << "Error in State::validState. The diffusion constant of the negative electrode is "
                  << st.Dn() << ", which is too high. The maximum is "
                  << s_ini.Dn() << ", the original effective diffusion constant.\n";
        throw 15;
    }

    // specific resistance of the electrodes (one value for both electrodes combined)
    bool rmin = st.r() <= 0;
    if (rmin)
    {
        std::cerr << "Error in State::validState. The specific resistance is "
                  << st.r() << ", which is too low, it must be strictly positive.\n";
        throw 15;
    }

    bool rmax = st.r() > 1000 * s_ini.r();
    if (rmax)
    {
        std::cerr << "Error in State::validState. The cell has degraded too much and the specific resistance is " << st.r()
                  << ", which is too high. The maximum is " << 1000 * s_ini.r() << ", 1000 times the original specific resistance.\n";
        throw 15;
    }

    // thickness of the plated litium layer
    bool delpl = st.delta_pl() < 0;
    if (delpl)
    {
        std::cerr << "Error in State::validState. The thickness of the plated lithium is " << st.delta_pl()
                  << ", which is too low. Only strictly positive values are allowed.\n";
        throw 15;

        // 0 is allowed
        // a negative value might lead to a further decrease in thickness, so it will keep getting more and more negative
    }

    // throw an error if one of the states was invalid
}