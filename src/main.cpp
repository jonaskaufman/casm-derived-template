#include <casm/crystallography/Structure.hh>
#include <casm/clex/ConfigMapping.hh>
#include <casm/clex/PrimClex.hh>
#include <casm/CASM_global_definitions.hh>
#include <casm/casm_io/VaspIO.hh>
#include <iostream>

CASM::Structure get_primitive(const CASM::Structure& input)
{
    CASM::Structure true_prim;
    input.is_primitive(true_prim);
    return true_prim;
}

double structure_score(const CASM::Structure& map_reference_struc,
                       const CASM::Structure& mappable_struc,
                       double weight = 0.0,
                       double vol_tol = 0.5,
                       double tol = CASM::TOL)
{
    // get prim and make PrimClex
    auto ref_prim=get_primitive(map_reference_struc);	
    CASM::Log log(std::cout, 0);
    CASM::Logging logging(log);
    CASM::PrimClex pclex(ref_prim, logging);
   
    // map it 
    CASM::Structure mappable_copy(mappable_struc);  // can't be const, make copy
    int options = 2;  // robust mapping
    CASM::ConfigMapper configmapper(pclex, weight, vol_tol, options, tol);
    CASM::jsonParser out;           // mapping output
    std::string name;               // not used
    std::vector<CASM::Index> best;  // not used
    Eigen::Matrix3d cart_op;        // not used
    bool update = false;
    configmapper.import_structure_occupation(mappable_copy, name, out, best, cart_op, update);
    double basis = out["best_mapping"]["basis_deformation"].get<double>();
    double lattice = out["best_mapping"]["lattice_deformation"].get<double>();
    return weight * lattice + (1 - weight) * basis;

}

int main()
{
    /**
     * Let's find the mapping score between two POSCARs
     */  

    CASM::Structure ref_struc("POSCAR_a");
    CASM::Structure map_struc("CONTCAR_a");
    std::cout << "Basis score: " << structure_score(ref_struc, map_struc, 0.0) << std::endl;
    std::cout << "Lattice score: " << structure_score(ref_struc, map_struc, 1.0) << std::endl;
    return 0;
}

