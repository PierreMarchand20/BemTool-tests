#ifndef BEMTOOLEXT_MISC_GMSHCALL_HPP
#define BEMTOOLEXT_MISC_GMSHCALL_HPP

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
		circle << "lc = "+NbrToStr(lc)+";\n";
		circle << "Point(0) = { 0 , 0 , 0 , lc};\n";

		//// Droite
		circle << "Point(1) = { "+NbrToStr(R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Haut
		circle << "Point(2) = { "+NbrToStr(0)+" , "+NbrToStr(R)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Gauche
		circle << "Point(3) = { "+NbrToStr(-R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Bas
		circle << "Point(4) = { "+NbrToStr(0)+" , "+NbrToStr(-R)+" , "+NbrToStr(0)+" , lc}; \n";

		//// Droite
		circle << "Circle(1) = { "+NbrToStr(1)+" , "+NbrToStr(0)+" , "+NbrToStr(2)+"}; \n";
		//// Haut
		circle << "Circle(2) = { "+NbrToStr(2)+" , "+NbrToStr(0)+" , "+NbrToStr(3)+"}; \n";
		//// Gauche
		circle << "Circle(3) = { "+NbrToStr(3)+" , "+NbrToStr(0)+" , "+NbrToStr(4)+"}; \n";
		//// Bas
		circle << "Circle(4) = { "+NbrToStr(4)+" , "+NbrToStr(0)+" , "+NbrToStr(1)+"}; \n";


		circle << "Physical Line(0) ={ "+NbrToStr(1)+" , "+NbrToStr(2)+" , "+NbrToStr(3)+" , "+NbrToStr(4)+"};\n";


		circle.close();
	}
	else std::cout << "Unable to open file \n" <<std::endl;

	system(("gmsh -2 -v "+NbrToStr(verbose)+" "+mesh_name+".geo").c_str());
}

////=============================================================////
////===========================  Disc =========================////
////=============================================================////
inline void gmsh_disc(std::string mesh_name, Real R, Real lc, int verbose=0){
	std::ofstream circle((mesh_name+".geo").c_str());
	if (circle.is_open()) {
		circle << "lc = "+NbrToStr(lc)+";\n";
		circle << "Point(0) = { 0 , 0 , 0 , lc};\n";

		//// Droite
		circle << "Point(1) = { "+NbrToStr(R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Haut
		circle << "Point(2) = { "+NbrToStr(0)+" , "+NbrToStr(R)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Gauche
		circle << "Point(3) = { "+NbrToStr(-R)+" , "+NbrToStr(0)+" , "+NbrToStr(0)+" , lc}; \n";
		//// Bas
		circle << "Point(4) = { "+NbrToStr(0)+" , "+NbrToStr(-R)+" , "+NbrToStr(0)+" , lc}; \n";

		//// Droite
		circle << "Circle(1) = { "+NbrToStr(1)+" , "+NbrToStr(0)+" , "+NbrToStr(2)+"}; \n";
		//// Haut
		circle << "Circle(2) = { "+NbrToStr(2)+" , "+NbrToStr(0)+" , "+NbrToStr(3)+"}; \n";
		//// Gauche
		circle << "Circle(3) = { "+NbrToStr(3)+" , "+NbrToStr(0)+" , "+NbrToStr(4)+"}; \n";
		//// Bas
		circle << "Circle(4) = { "+NbrToStr(4)+" , "+NbrToStr(0)+" , "+NbrToStr(1)+"}; \n";


		circle << "Line Loop(0) ={ "+NbrToStr(1)+" , "+NbrToStr(2)+" , "+NbrToStr(3)+" , "+NbrToStr(4)+"};\n";
		circle << "Plane Surface(0) = {"+NbrToStr(0)+"};\n";
		circle << "Physical Surface(0) = {"+NbrToStr(0)+"};\n";

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
