#include <iostream>
#include <complex>
#include <vector>

#include <htool/lrmat/partialACA.hpp>
#include <.h>

using namespace std;
using namespace htool;


int main(){
	const int ndistance = 4;
	double distance[ndistance];
	distance[0] = 15; distance[1] = 20; distance[2] = 30; distance[3] = 40;
	SetNdofPerElt(1);
	bool test = 0;
	for(int idist=0; idist<ndistance; idist++)
	{
		cout << "Distance between the clusters: " << distance[idist] << endl;
		SetEpsilon(0.0001);

        // TODO load meshes and create dof

        // TODO Creates clusters (needed EVEN in a low rank matrix)
		Cluster t(p1,tab1); Cluster s(p2,tab2);

        // generator
        BIO_Generator<HE_SL_3D_P1xP1,P1_2D> A(dof,kappa);

		// ACA with fixed rank
		int reqrank_max = 10;
		partialACA<double> A_partialACA_fixed(Ir,Ic,reqrank_max);
		A_partialACA_fixed.build(A,t,p1,tab1,s,p2,tab2);
		std::vector<double> partialACA_fixed_errors;
		for (int k = 0 ; k < A_partialACA_fixed.rank_of()+1 ; k++){
			partialACA_fixed_errors.push_back(Frobenius_absolute_error(A_partialACA_fixed,A,k));
		}

		cout << "Partial ACA with fixed rank" << endl;
		// Test rank
		cout << "rank : "<<A_partialACA_fixed.rank_of() << endl;
		test = test || !(A_partialACA_fixed.rank_of()==reqrank_max);

		// Test Frobenius errors
		test = test || !(partialACA_fixed_errors.back()<1e-8);
		cout << "Errors with Frobenius norm : "<<partialACA_fixed_errors<<endl;

		// Test compression
		test = test || !(0.87<A_partialACA_fixed.compression() && A_partialACA_fixed.compression()<0.89);
		cout << "Compression rate : "<<A_partialACA_fixed.compression()<<endl;

		// Test mat vec prod
		std::vector<double> f(nc,1);
		double error=norm2(A*f-A_partialACA_fixed*f);
		test = test || !(error<1e-6);
		cout << "Errors on a mat vec prod : "<< error<<endl<<endl;

		// ACA automatic building
		partialACA<double> A_partialACA(Ir,Ic);
		A_partialACA.build(A,t,p1,tab1,s,p2,tab2);
		std::vector<double> partialACA_errors;
		for (int k = 0 ; k < A_partialACA.rank_of()+1 ; k++){
			partialACA_errors.push_back(Frobenius_absolute_error(A_partialACA,A,k));
		}

		cout << "Partial ACA" << endl;
		// Test Frobenius error
		test = test || !(partialACA_errors[A_partialACA.rank_of()]<GetEpsilon());
		cout << "Errors with Frobenius norm: "<<partialACA_errors<<endl;

		// Test compression rate
		test = test || !(0.93<A_partialACA.compression() && A_partialACA.compression()<0.96);
		cout << "Compression rate : "<<A_partialACA.compression()<<endl;

		// Test mat vec prod
		error = norm2(A*f-A_partialACA*f);
		test = test || !(error<GetEpsilon()*10);
		cout << "Errors on a mat vec prod : "<< error<<endl<<endl<<endl;
	}
	cout << "test : "<<test<<endl;
	return test;
}
