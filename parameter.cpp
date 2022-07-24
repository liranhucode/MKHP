
#include "parameter.h"

using namespace std;

// ======================================================================= //
//  parameter
// ======================================================================= //
s_params::s_params()
{
	reduct_ratio = (float)1.7;
	coarseTo = 100;
}
// ======================================================================= //


// ======================================================================= //
s_params::~s_params()
{


}

// ======================================================================= //

// ======================================================================= //
//parse command: 
//				hypergraph file
//				numbers of partition
//              user load imbalance factor
//              numbers of iteration
// ======================================================================= //
void s_params::parse_command(int argc, char *argv[])
{

	int i;
	int num_read;

	//determine whether the number of command line parameters is correct
	if (argc < 5) {
		cout << "Usage : cmetis hypergraph.hgr nparts ubfactor nruns" << endl;
		exit(1);
	}

	filename = argv[1];


	num_read = sscanf_s(argv[2], "%d", &nparts);
	if (num_read != 1) {
		cout << "Error: %s option requires an integer parameter.\n", argv[2];
		exit(1);
	}

	num_read = sscanf_s(argv[3], "%d", &ubfactor);
	if (num_read != 1) {
		cout << "Error: %s option requires an integer parameter.\n", argv[3];
		exit(1);
	}

	num_read = sscanf_s(argv[4], "%d", &nruns);
	if (num_read != 1) {
		cout << "Error: %s option requires an integer parameter.\n", argv[4];
		exit(1);
	}



}
// ============================================================================================== //