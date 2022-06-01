#include "cantera/thermo.h"
#include <iostream>
#include "cantera/IdealGasMix.h"
#include "cantera/zerodim.h"
#include "vector"
#include <sstream>

//static CommDom* ptr_class = NULL;

// Declare some functions and pointes
//std::vector<std::unique_ptr<Cantera::IdealGasMix>> gasint;
//std::vector<std::unique_ptr<Cantera::IdealGasConstPressureReactor>> React;
//std::vector<std::unique_ptr<Cantera::ReactorNet>> net;

static Cantera::IdealGasMix* ptr_gas = NULL;

static Cantera::ReactorNet* ptr_net = NULL;

void cantera_sum (double* cc, double* cp_mass_frac, int cp_nsp, double* W_k,double* cp_low, double* cp_high); 

extern "C"
{

  const char*
  string_f2c(const char* licName, int nameLen)
  {
    char* name;
    name = (char*)malloc(sizeof(char)*(nameLen));
    for(int i = 0; i<nameLen; i++) name[i] = licName[i];
    name[nameLen] = '\0';
    return name ;
  }

    // Gets the Nasa coeffcients of each species (Need to be run only once)
    void cantera_coeff_(const char* fmech, int* n_fmech, double* cc, double* W_k)
    {   int type, nsp;
        double temp[15]; 
        double minTemp, maxTemp, refPressure;

        Cantera::IdealGasMix gas(string_f2c(fmech,n_fmech[0]));
        nsp = gas.nSpecies();
        gas.getMolecularWeights(W_k);
        for(int hh=0; hh<nsp; hh++)  
            W_k[hh] = W_k[hh] / 1000.0; // Convert into SI units (kg/mol)

        Cantera::MultiSpeciesThermo& sp = gas.speciesThermo();
        for(int jj=0,kk=0; jj<nsp; jj++){ 
            sp.reportParams(jj, type, temp, minTemp, maxTemp, refPressure);
            for(int ii=0; ii<15; ii++,kk++){
                cc[kk] = temp[ii];
            }
        }
    }

    void cantera_divide_(double *num, double *den, double *res)
    {
       *res = *num/ *den;
    }

    void cantera_log_(double *in)
    {
       *in = log10(*in);
    }
    void cantera_alya_cp_(int* nsp, double* cc, double  *mass_frac, double* W_k , double *low_coeffs, double *high_coeffs)
    {
        cantera_sum (cc, mass_frac, nsp[0], W_k ,low_coeffs, high_coeffs);
    }

    void cantera_elemh_(double *mass_frac, double *H)
    {   int ha;
	ptr_gas[0].setMassFractions(&mass_frac[0]);
        ha = ptr_gas[0].elementIndex("H");
        *H = ptr_gas[0].elementalMassFraction(ha);
    }
    void cantera_elemc_(double *mass_frac, double *C)
    {   int ca;
	ptr_gas[0].setMassFractions(&mass_frac[0]);
        ca = ptr_gas[0].elementIndex("C");
        *C = ptr_gas[0].elementalMassFraction(ca);

    }
    void cantera_elemo_(double *mass_frac, double *O)
    {   int oa;
	ptr_gas[0].setMassFractions(&mass_frac[0]);
        oa = ptr_gas[0].elementIndex("O");
        *O = ptr_gas[0].elementalMassFraction(oa);
    }
    void cantera_reduced_(char* non_key)
    {   
        int size_red;
        std::string str1 ("END_REDUCTION"); 
        int reac_multi[ptr_gas[0].nReactions()];
        for (int j=0; j<ptr_gas[0].nReactions(); j++) reac_multi[j] = 1;
        // Delimit "non_key" string to get an array of non key species
        std::string s;
        std::stringstream ss;
        ss << non_key;
        ss >> s;
        std::vector<std::string>   result;
        std::stringstream  data(s);
        std::string line;
        while(std::getline(data,line,','))
        {
            result.push_back(line);
        }
        size_red = result.size();
        // Get index of reactions including the non key species
        int spec_index[size_red];
        if (str1.compare(result[0]) != 0) {
            for (int j=0; j<size_red; j++) spec_index[j] = ptr_gas[0].speciesIndex(result[j]);
            for(int i=0; i<size_red; i++){
                for(int j=0; j < ptr_gas[0].nReactions(); j++){
                    if ( (ptr_gas[0].reactantStoichCoeff(spec_index[i], j) > 0.001 || ptr_gas[0].productStoichCoeff(spec_index[i], j) > 0.001) )
                        reac_multi[j] = 0;
                }
            }
            for(int j=0; j < ptr_gas[0].nReactions(); j++) ptr_gas[0].setMultiplier(j, reac_multi[j]);
        }
        result.clear();
    }

    void cantera_prog_(char* prog, int *prog_index)
    {   
        int size_red;
        // Delimit "non_key" string to get an array of non key species
        std::string s;
        std::stringstream ss;
        ss << prog;
        ss >> s;
        std::vector<std::string>   result;
        std::stringstream  data(s);
        std::string line;
        while(std::getline(data,line,','))
        {
            result.push_back(line);
        }
        size_red = result.size();
        for (int j=0; j<size_red; j++) prog_index[j]  = ptr_gas[0].speciesIndex(result[j]);
        result.clear();
    }

    // Integration Source terms
    void cantera_initialization_(int *n_size, const char* fmech, int* n_fmech)
    {
        ptr_gas = new Cantera::IdealGasMix(string_f2c(fmech,n_fmech[0]));
    }

    void cantera_trim_(const char* fmech, int* n_fmech, int* key_index)
    {
        int size_key;   
        std::string s;
        std::stringstream ss;
        ss << string_f2c( fmech,n_fmech[0] );
        ss >> s;
        std::vector<std::string>   result;
        std::stringstream  data(s);
        std::string line;
        while(std::getline(data,line,','))
        {
            result.push_back(line);
        }
        size_key = result.size();
        for(int j=0; j < size_key; j++) key_index[j] = ptr_gas[0].speciesIndex(result[j]);
    }

    void cantera_dac_integrate_(double *temp, double *p,  double *mass_frac, int* reac_multi, double* t)
    {
        for(int j=0; j < ptr_gas[0].nReactions(); j++) ptr_gas[0].setMultiplier(j, reac_multi[j]);
        ptr_gas[0].setState_TPY(*temp, *p, &mass_frac[0]);
        Cantera::IdealGasConstPressureReactor r;
        r.insert(ptr_gas[0]);
        Cantera::ReactorNet sim;
        sim.addReactor(r);
        // Integrate Source Terms
        sim.advance(*t);
        *temp =  r.temperature();
        *p =  r.pressure();
        ptr_gas[0].getMassFractions(&mass_frac[0]);
    }

    void cantera_integrate_(double *temp, double *p,  double *mass_frac, double *t)
    { 
        ptr_gas[0].setState_TPY(*temp, *p, &mass_frac[0]);
        Cantera::IdealGasConstPressureReactor r;
        r.insert(ptr_gas[0]);
        Cantera::ReactorNet sim;
        sim.addReactor(r);
        // Integrate Source Terms
        sim.advance(*t);
        *temp =  r.temperature();
        *p =  r.pressure();
        ptr_gas[0].getMassFractions(&mass_frac[0]);
    }

    void cantera_equi_integrate_(double *temp, double *p,  double *mass_frac, double *mass_frac_eq, double *t, int* all_spec_index)
    { 
        for(int j=0; j < ptr_gas[0].nSpecies(); j++){
            if(all_spec_index[j] == 0) mass_frac[j] = mass_frac_eq[j];
        }
        ptr_gas[0].setState_TPY(*temp, *p, &mass_frac[0]);
        Cantera::IdealGasConstPressureReactor r;
        r.insert(ptr_gas[0]);
        Cantera::ReactorNet sim;
        sim.addReactor(r);
        // Integrate Source Terms
        sim.advance(*t);
        *temp =  r.temperature();
        *p =  r.pressure();
        ptr_gas[0].getMassFractions(&mass_frac[0]);
    }

    void cantera_equilibrate_(double *temp, double *p,  double *mass_frac,char* non_key, int* all_spec_index)
    {   
        int size_red;
        int reac_multi[ptr_gas[0].nReactions()]; 
        for (int j=0; j<ptr_gas[0].nReactions(); j++) reac_multi[j] = 1;
        // Delimit "non_key" string to get an array of non key species
        std::string s;
        std::stringstream ss;
        ss << non_key;
        ss >> s;
        std::vector<std::string>   result;
        std::stringstream  data(s);
        std::string line;
        while(std::getline(data,line,','))
        {
            result.push_back(line);
        }
        size_red = result.size();
        // Get index of reactions including the non key species
        int spec_index[size_red];
        for (int j=0; j<size_red; j++) spec_index[j] = ptr_gas[0].speciesIndex(result[j]);

        for (int j=0; j<ptr_gas[0].nSpecies(); j++) {
            for (int i=0; i<size_red; i++) {
                if (j == spec_index[i]) all_spec_index[j] = 0;
            }
        }
        ptr_gas[0].setState_TPY(*temp, *p, &mass_frac[0]);
        ptr_gas[0].equilibrate("HP");
        ptr_gas[0].getMassFractions(&mass_frac[0]); 
        result.clear();
    }

    void cantera_partial_molar_enthalpies_(int *nspec, double *temp, double *p, double *mass_frac, double *partial_h)
    {
        ptr_gas[0].setState_TPY(*temp, *p, &mass_frac[0]);
        ptr_gas[0].getPartialMolarEnthalpies(&partial_h[0]);
        for(int ii=0; ii < *nspec; ii++)
            partial_h[ii] = partial_h[ii] / 1000.0;  // Convert into SI units (J/mol)
    }

    void cantera_get_enthalpy_from_ty_(double *p, double *temp, double *mass_frac, double *enthalpy_spec_RT)
    {
        ptr_gas[0].setState_TPY(*temp, *p, &mass_frac[0]);
        ptr_gas[0].getEnthalpy_RT(&enthalpy_spec_RT[0]);
    }

    void cantera_get_temperature_from_hy_(double *p, double *enthalpy, double *mass_frac, double *temp)
    {
        ptr_gas[0].setMassFractions(&mass_frac[0]);
        ptr_gas[0].setState_HP(*enthalpy, *p);
        *temp = ptr_gas[0].temperature();
    }

    void cantera_equilibrium_from_hy_(double *p, double *enthalpy, double *mass_frac, double *temp_eq)
    {
        ptr_gas[0].setMassFractions(&mass_frac[0]);
        ptr_gas[0].setState_HP(*enthalpy, *p);
        ptr_gas[0].equilibrate("HP");
        ptr_gas[0].getMassFractions(&mass_frac[0]);
        *temp_eq = ptr_gas[0].temperature();
    }

} // extern "C"

// sums up the matrix rows to give the mixture coeffceints 
void cantera_sum (double* cp_cc,double* cp_mass_frac, int cp_nsp, double* W_k ,double* cp_low, double* cp_high)
{
    double sum_coeffs[15];
    std::vector< std::vector<double> > xxx(cp_nsp);
    for(int ii=0; ii<cp_nsp; ii++) xxx[ii] =  std::vector<double>(15,0);

    for(int ii=0, gg=0; ii<cp_nsp; ii++) {	 
        for(int jj=0; jj<15; jj++, gg++) xxx[ii][jj] = (cp_cc[gg]*Cantera::GasConstant*cp_mass_frac[ii]/W_k[ii]) / 1000.0;
    }
    for(int ii=0; ii<15; ii++) 
    {sum_coeffs[ii] = 0;
        for(int jj=0; jj<cp_nsp; jj++) {
            sum_coeffs[ii] = sum_coeffs[ii] + xxx[jj][ii];
        }
    }
    for(int ii=0; ii<6; ii++) cp_low[ii]  = sum_coeffs[8+ii]; 
    for(int ii=0; ii<6; ii++) cp_high[ii] = sum_coeffs[1+ii];
}
