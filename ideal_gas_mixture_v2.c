/*
 * =============================================================================
 *  File:        real_ideal.c
 *  Original ID: @(#)real_ideal.c    1.10
 *  Original Copyright © 1988–1998 ANSYS, Inc.
 *  All Rights Reserved.
 *
 *  This file contains unpublished proprietary source code of ANSYS, Inc.
 *  It is protected by U.S. copyright law as an unpublished work and is
 *  furnished under a written license agreement. It is considered confidential
 *  and may not be used, copied, modified, or disclosed except in accordance
 *  with the terms of that agreement.
 * =============================================================================
 *
 *  Forked and extended version by Arash Jalil-Khabbazi
 *  Copyright © 2023 Arash Jalil-Khabbazi
 *
 *  This fork builds upon the original ANSYS framework to implement and
 *  demonstrate property calculations for ideal gas mixtures within the
 *  ANSYS Fluent environment. Intended for academic and illustrative use only.
 *
 * =============================================================================
 *
 *  ⚠ WARNING:
 *  Including "udf.h" is solely for accessing definitions of ANSYS Fluent
 *  constructs (e.g., Domain). You must NOT reference internal Fluent global
 *  variables or link against Fluent binaries (e.g., fl551.exe), as this
 *  creates version-dependent behavior and violates Fluent’s modular UDF design.
 * =============================================================================
 */


/* ALL VARIABLES ARE IN SI UNITS.
 *    
 * 	pressure [Pa]
 * 	temperature [K]
 * 	molar_weight [kg/kmol]
 *  specific_volume [m3/kg]
 *  molar_volume [m3/kmol]
 *  universal_gas_constant (molar_gas_constant) [J/kmol-K]
 *  mass_specific_gas_constant (J/kg-K)
 *  enthalpy [J/kg]
 *  entropy [J/kg-K]
 *  specific_heat [J/kg-K]
 *
 */ 

#include "udf.h"
#include "stdio.h"
#include "ctype.h"
#include "stdarg.h"

#if RP_DOUBLE
#define SMLL 1.e-20
#else
#define SMLL 1.e-10
#endif

#define NSPECIE_NAME 80 
#define RGASU UNIVERSAL_GAS_CONSTANT /* 8314.34 SI units: J/kmol/K */
#define PI  3.141592654

/* Here input the number of species in the mixture */
/* THIS IS A USER INPUT */
#define n_specs 2

/* OPTIONAL REFERENCE (OFFSET) VALUES FOR ENTHALPY AND ENTROPY  */
#define h_ref 0.0
#define s_ref 0.0

static int (*usersMessage)(const char *,...);
static void (*usersError)(const char *,...);

/* REFERENCE STATE */
#define P_ref 101325.0
#define T_ref 298.15
//static double P_ref, T_ref;
static char  gas[n_specs][NSPECIE_NAME];

/* static property parameters */
static double cp[5][n_specs]; /* specific heat polynomial coefficients */
static double mw[n_specs];  /* molecular weights */
static double tcrit[n_specs]; /* critical temperature */
static double pcrit[n_specs]; /* critical pressure */
static double vcrit[n_specs]; /* critical specific volume */

/* Static variables associated with ideal gas model */
static double rgas[n_specs], cp_inT_ref[n_specs], h_ideal_ref[n_specs];

void Mw();
void Cp_Parameters(); 
void Tcrit();
void Pcrit();
void Vcrit();

double cp_ideal_gas(double temp, int i);
double viscosity(double temp, int i);
double thermal_conductivity(double temp, int i);

DEFINE_ON_DEMAND(I_do_nothing)
{
    /*
    This is a dummy function
    must be included to allow for the use of the
    ANSYS FLUENT UDF compilation utility
    */
}

/*******************************************************************/
/* Mixture Functions                                               */ 
/* These are the only functions called from ANSYS FLUENT Code      */ 
/*******************************************************************/

void MIXTURE_Setup(Domain *domain, cxboolean vapor_phase, char *specielist,
                    int (*messagefunc)(const char *format,...),
                    void (*errorfunc)(const char *format,...))
{
    /* This function will be called from ANSYS FLUENT after the
    UDF library has been loaded.
    User must enter the number of species in the mixture
    and the name of the individual species.
    */
	
	int i, j;
	usersMessage = messagefunc;
    usersError = errorfunc;

	/*
    P_ref = ABS_P(RP_Get_Real("reference-pressure"), op_pres);
    T_ref = 298.15;

    Message0("\n MIXTURE_Setup: ideal gas equation of State \n");
    Message0("\n MIXTURE_Setup: reference-temperature is %f \n", T_ref);

    if (P_ref == 0.0){
        Message0("\n MIXTURE_Setup: reference-pressure was not set by user \n");
        Message0("\n MIXTURE_Setup: setting reference-pressure to 101325 Pa \n");
        P_ref = 101325.0;
    }
	*/
	
    /*====================================================*/
    /*=========  User Input Section ======================*/
    /*====================================================*/
    /*
    Define Species name.
    DO NOT use space for naming species
    */
    (void)strcpy(gas[0], "H2"); 
    (void)strcpy(gas[1], "CH4");
 
    /*====================================================*/
    /*=========  End Of User Input Section ===============*/
    /*====================================================*/
	
    Message0("\n MIXTURE_Setup: RealGas mixture initialization \n");
    Message0("\n MIXTURE_Setup: Number of Species = %d \n", n_specs);
    for (i=0; i<n_specs; ++i)
    {   
        Message0("\n MIXTURE_Setup: Species[%d] = %s \n", i, gas[i]);
    }

    /*
    concatenate species name into one string and send back to fluent
    */
    strcat(specielist, gas[0]);
    for (i=1; i<n_specs; ++i)
    {
        strcat(specielist, " ");
        strcat(specielist, gas[i]);
    }

    /* initialize */
    Mw();
    Cp_Parameters(); 
    Tcrit();
    Pcrit();
    Vcrit();

    for (i=0; i<n_specs; ++i)
    {
        rgas[i] = RGASU/mw[i];
        h_ideal_ref[i] = T_ref*(cp[0][i] + T_ref*(0.5*cp[1][i] + T_ref*(1./3.*cp[2][i] + T_ref*(0.25*cp[3][i] + T_ref*0.2*cp[4][i]))));
        cp_inT_ref[i] = cp[0][i]*log(T_ref)+T_ref*(cp[1][i]+T_ref*(0.5*cp[2][i]+T_ref*(1./3.*cp[3][i]+0.25*cp[4][i]*T_ref))); 
    }
    
}

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_mw                                       */
/* Returns mixture molar mass                                 */
/*------------------------------------------------------------*/

double MIXTURE_mw(double yi[])
{
    double sum = 0.0;
    int i;

    for (i=0; i<n_specs; ++i)
        sum += yi[i]/mw[i];
 
    return 1.0/sum; /* (kg/kmol) */
}

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_rgas                                     */
/* Returns mass-specific mixture gas constant                 */
/*------------------------------------------------------------*/

double MIXTURE_rgas(double yi[])
{
    double rgas_m = 0.0;
    int i;
    for (i=0; i<n_specs; ++i)
        rgas_m += yi[i]*rgas[i];
    
    return rgas_m; /* (J/kg/K) */
}

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_density                                  */
/*      Returns density given T and P                         */
/*------------------------------------------------------------*/

double MIXTURE_density(cell_t cell, Thread *thread, cxboolean vapor_phase, double temp, double P, double yi[])
{
    return P/(MIXTURE_rgas(yi)*temp); /* (kg/m^3) */
}

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_specific_heat                            */
/*      Returns specific heat given T and rho                 */
/*------------------------------------------------------------*/

double MIXTURE_specific_heat(cell_t cell, Thread *thread, double temp, double density, double P, double yi[])
{
    double cp = 0.0;
	int i;
	
    for (i=0; i<n_specs; ++i)
        cp += yi[i]*(cp[0][i] + temp*(cp[1][i] + temp*(cp[2][i] + temp*(cp[3][i] + temp*cp[4][i]))));

    return cp; /* (J/kg/K) */
}

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_enthalpy                                 */
/*      Returns specific enthalpy given T and rho             */
/*------------------------------------------------------------*/

double MIXTURE_enthalpy(cell_t cell, Thread *thread, double temp, double density, double P, double yi[])
{		  
    double h = 0.0;
	int i;
		
    for (i=0; i<n_specs; ++i)
        h += yi[i]*(temp*(cp[0][i] + temp*(0.5*cp[1][i] + temp*(1./3.*cp[2][i] + temp*(0.25*cp[3][i] + temp*0.2*cp[4][i]))))
        - h_ideal_ref[i]);
    
    return h + h_ref; /* (J/kg) */

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_cp_integral                              */
/*      Returns entropy's Cp_integral part given T and rho    */
/*------------------------------------------------------------*/

double MIXTURE_cp_integral(double temp, double yi[])
{
    double cp_int = 0.0;
	int i;

    for (i=0; i<n_specs; ++i)
        cp_int += yi[i]*(cp[0][i]*log(temp)+temp*(cp[1][i]+temp*(0.5*cp[2][i]+temp*(1./3.*cp[3][i]+0.25*cp[4][i]*temp)))
        - cp_inT_ref[i]);

    return cp_int; /* (J/kg-K) */
}

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_entropy                                  */
/*      Returns entropy given T and rho                       */
/*------------------------------------------------------------*/

double MIXTURE_entropy(cell_t cell, Thread *thread, double temp, double density, double P, double yi[])
{
    return (s_ref + MIXTURE_cp_integral(temp, yi) + MIXTURE_rgas(yi)*log(P_ref/P)); /* (J/kg/K) */
}

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_speed_of_sound                           */
/*      Returns s.o.s given T and rho                         */
/*------------------------------------------------------------*/

double MIXTURE_speed_of_sound(cell_t cell, Thread *thread, double temp, double density, double P, double yi[])
{
    double cp = MIXTURE_specific_heat(cell, thread, temp, density, P, yi);
    double rgas = MIXTURE_rgas(yi);
    
    return sqrt(rgas*temp*cp/(cp-rgas)); /* m/s */
}

/*------------------------------------------------------------------------*/
/* FUNCTION: MIXTURE_rho_t                                                */
/* derivative of denisty w.r.t. temperature at constant pressure (drdT|p) */
/*------------------------------------------------------------------------*/

double MIXTURE_rho_t(cell_t cell, Thread *thread, double temp, double density, double P, double yi[])
{
    return -density/temp; /* (kg/m^3/K) */
}

/*------------------------------------------------------------------------*/
/* FUNCTION: MIXTURE_rho_p                                                */
/* derivative of denisty w.r.t. pressure at constant temperature (drdp|T) */
/*------------------------------------------------------------------------*/

double MIXTURE_rho_p(cell_t cell, Thread *thread, double temp, double density, double P, double yi[])
{
    double rgas = MIXTURE_rgas(yi);
    return 1./(rgas*temp); /* (kg/m^3/Pa) */
}

/*-------------------------------------------------------------------------*/
/* FUNCTION: MIXTURE_enthalpy_t                                            */
/* derivative of enthalpy w.r.t. temperature at constant pressure (dhdT|p) */
/*-------------------------------------------------------------------------*/

double MIXTURE_enthalpy_t(cell_t cell, Thread *thread, double Temp, double density, double P, double yi[])
{
    return MIXTURE_specific_heat(cell, thread, Temp, density, P, yi); /* J/(kg.K) */
}

/*-------------------------------------------------------------------------*/
/* FUNCTION: MIXTURE_enthalpy_p                                            */
/* derivative of enthalpy w.r.t. pressure at constant temperature (dhdp|T) */
/*-------------------------------------------------------------------------*/

double MIXTURE_enthalpy_p(cell_t cell, Thread *thread, double Temp, double density, double P, double yi[])
{
    /* general form dh/dp|T = (1/rho)*[ 1 + (T/rho)*drho/dT|p] */
    /* but for ideal gas dh/dp = 0 */
    return 0.0; /* J/(kg.Pa) */
}

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_viscosity                                */
/*------------------------------------------------------------*/

double MIXTURE_viscosity(cell_t cell, Thread *thread, double temp, double density, double P, double yi[])
{
    double mu = 0.;
	int i;
	
	for (i=0; i<n_specs; ++i)
		mu += yi[i]*viscosity(temp, i);
	
	return mu; /* kg/(m.s) */
}

/*------------------------------------------------------------*/
/* FUNCTION: MIXTURE_thermal_conductivity                     */
/*------------------------------------------------------------*/

double MIXTURE_thermal_conductivity(cell_t cell, Thread *thread, double temp, double density, double P, double yi[])
{
    double kt = 0.;
	int i;
	
	for (i=0; i<n_specs; ++i)
        kt += yi[i]*thermal_conductivity(temp, i);
	
	return kt; /* W/(m.K) */
}

/*******************************************************************/
/* Species Property Definitions                                    */ 
/*******************************************************************/

void Mw() /* molecular weight */
{ /* kg/kmol */
    mw[0] = 2.01594; /*H2*/
    mw[1] = 16.04303; /*CH4*/
}

void Pcrit() /* critical pressure */
{ /* Pa */
    pcrit[0] = 1.293e6; /*H2*/
    pcrit[1] = 4.599e6; /*CH4*/
}

void Tcrit() /* critical temperature */
{ /* K */
    tcrit[0] = 32.98; /*H2*/
    tcrit[1] = 190.56; /*CH4*/
}

void Vcrit() /* critical specific volume */
{ /* m3/kg */
    vcrit[0] = 0.031846; /*H2*/
    vcrit[1] = 0.006146; /*CH4*/
}

void Cp_Parameters() /* coefficients of specific heat polynomials */
{ /* J/kg/K */
    cp[0][0] = 13602.45 ; /*H2*/
    cp[1][0] = 3.402317; 
    cp[2][0] = -0.003358423;
    cp[3][0] = -3.907953e-07;
    cp[4][0] = 1.705345e-09;

    cp[0][1] = 403.5847 ; /*CH4*/
    cp[1][1] = 9.057335; 
    cp[2][1] = -0.01442509;
    cp[3][1] = 1.580519e-05;
    cp[4][1] = -6.343051e-09 ;
}

/*------------------------------------------------------------*/
/* FUNCTION: cp_ideal_gas                                     */
/*             Returns ideal gas specific heat given T        */
/*------------------------------------------------------------*/

double cp_ideal_gas(double temp, int i)
{	
	double cpi = (cp[0][i]+temp*(cp[1][i]+temp*(cp[2][i]+temp*(cp[3][i] +temp*cp[4][i]))));
	
	return cpi;
}

/*******************************************************************/
/* FUNCTION: viscosity                                             */
/*******************************************************************/

double viscosity(double temp, int i)
{
	double mu, tr, tc, pcatm;
	
	tr = temp/tcrit[i];
	tc = tcrit[i];
	pcatm = pcrit[i]/101325.;
	mu = 6.3e-7 * sqrt(mw[i])*pow(pcatm,0.6666)/pow(tc,0.16666) *(pow(tr,1.5)/(tr+0.8));
	
	return mu; /* (kg/m/s) */
}

/*******************************************************************/
/* FUNCTION: thermal_conductivity                           */
/*******************************************************************/

double thermal_conductivity(double temp, int i)
{
	double cp, mu;
	
	cp = cp_ideal_gas(temp, i);
	mu = viscosity(temp, i);

	return (cp + 1.25*rgas[i])*mu; /* W/m/K */
}

/*******************************************************************/
/* Mixture Functions Structure                                     */ 
/*******************************************************************/

UDF_EXPORT RGAS_Functions RealGasFunctionList = 
{
    MIXTURE_Setup,                    /* initialize            */
    MIXTURE_density,                  /* density               */
    MIXTURE_enthalpy,                 /* sensible enthalpy     */
    MIXTURE_entropy,                  /* entropy               */
    MIXTURE_specific_heat,            /* specific_heat         */
    MIXTURE_mw,                       /* molecular_weight      */
    MIXTURE_speed_of_sound,           /* speed_of_sound        */
    MIXTURE_viscosity,                /* viscosity             */
    MIXTURE_thermal_conductivity,     /* thermal_conductivity  */
    MIXTURE_rho_t,                    /* drho/dT |const p      */
    MIXTURE_rho_p,                    /* drho/dp |const T      */
    MIXTURE_enthalpy_t,               /* dh/dT |const p        */
    MIXTURE_enthalpy_p,               /* dh/dp |const T        */
};
