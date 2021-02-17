#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <iostream> 
#include<cmath>
#include<stdio.h>
#include<iomanip>
// #include <boost/filesystem.hpp>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
//#include "skin_surface_writer.h"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::FT                                  FT;
typedef Skin_surface_3::Bare_point                          Bare_point;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef CGAL::Polyhedron_3<K,
  CGAL::Skin_surface_polyhedral_items_3<Skin_surface_3> >   Polyhedron;



using namespace std;

std::vector<float> read_csv(std::string filename){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>

    // Create a vector of <string, int vector> pairs to store the result
    std::vector<std::pair<std::string, std::vector<float>>> result;
    std::vector<float> resultf;

    //cout << "input-file_csv= " <<   filename  << "  \n"; 
    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname;
    float val;

    // Read the column names
    if(myFile.good())
    {
        // Extract the first line in the file
        std::getline(myFile, line);

        // Create a stringstream from line
        std::stringstream ss(line);

        // Extract each column name
        while(std::getline(ss, colname, ',')){
            
            // Initialize and add <colname, int vector> pairs to result
            result.push_back({colname, std::vector<float> {}});
            //std::cout << "header result: "  << colname << "\n"; 
        }
    }

    

    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);
        
        // Keep track of the current column index
        int colIdx = 0;
        
        // Extract each integer
        while(ss >> val){
            
            // Add the current integer to the 'colIdx' column's values vector
            result.at(colIdx).second.push_back(val);
            resultf.push_back(val);
            //std::cout << "col result: "  << colIdx << val << "\n"; 
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
            
            // Increment the column index
            colIdx++;
        }
    }

    // Close file 
    myFile.close();


    return resultf;
}


std::string format_cell_label(float cell_label) {
    ostringstream out;
    out << std::internal << std::setfill('0') << std::setw(3) << cell_label;
    return out.str();
}

void print_row_data(std::vector<float> v_row){
    cout << "cell_label " << v_row[0] << " : \n"; 
    cout << "x: " << v_row[1] << "\n"; 
    cout << "y: " << v_row[2] << "\n"; 
    cout << "z: " << v_row[3] << "\n"; 
    cout << "radius: " << v_row[4] << "\n"; 
    cout << "weight: " << pow(v_row[4],2) << "\n"; //A weighted point p=(p,wp)∈R3×R corresponds to a ball with center p and radius= sqrt(wp)


}

void mesh_cell(std::list<Weighted_point> l_wpoints, float cell_label, std::string output_folder, float shrinkf){
    //cout << "CGAL meshing  " << cell_label << " \n"; 
    FT shrinkfactor = shrinkf;
    Polyhedron p;

    Skin_surface_3 skin_surface(l_wpoints.begin(), l_wpoints.end(), shrinkfactor);
    CGAL::mesh_skin_surface_3(skin_surface, p);
    CGAL::subdivide_skin_surface_mesh_3(skin_surface, p);

    // std::ostringstream ss;
    // ss << cell_label;
    // std::string s_label(ss.str());

    std::string s_label(format_cell_label(cell_label));
    //cout << "cell_label is " << s_label << " : \n"; 
    std::string file_name = output_folder + "/" + "cell_" + s_label + ".off";
    std::ofstream out(file_name);
    out << p;
}

void mesh_spheres(std::vector<float> v_spheres, std::string output_folder, float shrinkfactor){
    //std::cout << "inside mesh_spheres";
    int colIdx = 0;
    float cell_label_prev = -1;
    std::vector<float> v_row;
    std::list<Weighted_point> l_wpoints;

    for (std::vector<float>::const_iterator i = v_spheres.begin(); i != v_spheres.end(); ++i) {
        //std::cout << *i << ' ';
        v_row.push_back(*i);
        if (colIdx == 0 && cell_label_prev != -1){
            if (cell_label_prev !=  v_row[0]){
                mesh_cell(l_wpoints,cell_label_prev,output_folder,shrinkfactor);
                l_wpoints.clear();
            }
        }

        if (colIdx == 4){
            //print_row_data(v_row);
            l_wpoints.push_front(Weighted_point(Bare_point( v_row[1],v_row[2],v_row[3]), pow(v_row[4],2))); 
            v_row.clear();
            colIdx = 0;
            cell_label_prev = v_row[0];
        }
        else {
            colIdx++;
        }
        

    }

    if (l_wpoints.size()>0){
        mesh_cell(l_wpoints,cell_label_prev,output_folder,shrinkfactor);

    }
}



int main(int argc, char *argv[]) {

    string input_file_csv, output_folder;
    float shrinkfactor;

    if(argc>=2) {
        input_file_csv = argv[1];
        output_folder = argv[2];
        shrinkfactor = (float)strtod(argv[3],NULL);

    }
    else{
        input_file_csv = "spheres_for_meshing_CGAL.csv";
        output_folder = "." ;
        shrinkfactor = 0.3;
    }


    // const char* path = output_folder.c_str();
    // boost::filesystem::path dir(path);
    // if(boost::filesystem::create_directory(dir))
    // {
    //     std::cerr<< "Directory Created: "<< output_folder <<std::endl;
    // }

    cout << "input-file_csv= " <<   input_file_csv <<  "output folder= " << output_folder << "  \n"; 
    std::vector<float> v_spheres = read_csv(input_file_csv);
    mesh_spheres(v_spheres,output_folder,shrinkfactor);

    return 0;
}
