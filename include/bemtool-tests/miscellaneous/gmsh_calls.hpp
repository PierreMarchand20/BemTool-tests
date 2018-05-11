#ifndef BEMTOOLTESTS_MISC_GMSHCALL_HPP
#define BEMTOOLTESTS_MISC_GMSHCALL_HPP

#include <string>
#include <fstream>
#include <bemtool/calculus/calculus.hpp>
#include <bemtool/miscellaneous/misc.hpp>

namespace bemtool{



////=============================================================////
////===========================  Circle =========================////
////=============================================================////
inline void gmsh_circle(std::string mesh_name, Real R, Real lc, int verbose=0){
	std::ofstream circle((mesh_name+".geo").c_str());
    if (circle.is_open()) {
        circle << "SetFactory(\"OpenCASCADE\");"<<std::endl;
        circle << "h = "<<NbrToStr(lc)<<";"<<std::endl;
        circle << "Mesh.CharacteristicLengthMin = h;"<<std::endl;
        circle << "Mesh.CharacteristicLengthMax = h;"<<std::endl;
        circle << "Circle(1) = {0,0,0,"<< NbrToStr(R)<<"};";
        circle.close();
    }
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -1 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////===========================  Disc =========================////
////=============================================================////
inline void gmsh_disc(std::string mesh_name, Real R, Real lc, int verbose=0){
	std::ofstream circle((mesh_name+".geo").c_str());
    if (circle.is_open()) {
        circle << "SetFactory(\"OpenCASCADE\");"<<std::endl;
        circle << "h = "<<NbrToStr(lc)<<";"<<std::endl;
        circle << "Mesh.CharacteristicLengthMin = h;"<<std::endl;
        circle << "Mesh.CharacteristicLengthMax = h;"<<std::endl;
        circle << "Disk(1) = {0,0,0,"<< NbrToStr(R)<<"};";
        circle.close();
    }
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -2 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////===========================  Sphere =========================////
////=============================================================////
inline void gmsh_sphere(std::string mesh_name, Real R, Real lc, int verbose=0){
	std::ofstream sphere((mesh_name+".geo").c_str());
	if (sphere.is_open()) {
		sphere << "SetFactory(\"OpenCASCADE\");"<<std::endl;
		sphere << "h = "<<NbrToStr(lc)<<";"<<std::endl;
		sphere << "Mesh.CharacteristicLengthMin = h;"<<std::endl;
		sphere << "Mesh.CharacteristicLengthMax = h;"<<std::endl;
		sphere << "Sphere(1) = {0,0,0,"<< NbrToStr(R)<<"};";
		sphere.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -2 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////===========================  Ball ===========================////
////=============================================================////
inline void gmsh_ball(std::string mesh_name, Real R, Real lc, int verbose=0){
	std::ofstream ball((mesh_name+".geo").c_str());
	if (ball.is_open()) {
		ball << "SetFactory(\"OpenCASCADE\");"<<std::endl;
		ball << "h = "<<NbrToStr(lc)<<";"<<std::endl;
		ball << "Mesh.CharacteristicLengthMin = h;"<<std::endl;
		ball << "Mesh.CharacteristicLengthMax = h;"<<std::endl;
		ball << "Sphere(1) = {0,0,0,"<< NbrToStr(R)<<"};";
		ball.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -3 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////==========================  Segment =========================////
////=============================================================////
inline void gmsh_seg(std::string mesh_name, R3 A, R3 B, int verbose=0){
	std::ofstream segment((mesh_name+".geo").c_str());
	if (segment.is_open()) {
		Real lc =0.1;
		segment << "lc = "+NbrToStr(lc)+";\n";

		//// Droite
		segment << "Point(0) = { "+NbrToStr(A[0])+" , "+NbrToStr(A[1])+" , "+NbrToStr(A[2])+" , lc}; \n";
		//// Haut
		segment << "Point(1) = { "+NbrToStr(B[0])+" , "+NbrToStr(B[1])+" , "+NbrToStr(B[2])+" , lc}; \n";

		segment << "Line(2) =  {0,1}; \n";

		segment << "Physical Line(3) ={3};\n";


		segment.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -2 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////===========================  Clean ==========================////
////=============================================================////
inline void gmsh_clean(std::string mesh_name){
	system(("rm "+mesh_name+".geo").c_str());
	system(("rm "+mesh_name+".msh").c_str());
}

}
#endif
