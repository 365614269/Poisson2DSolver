#include "include/MeshRect.h"
#include <fstream>

using namespace std;
using namespace rapidxml;
using Eigen::VectorXd;
using Eigen::MatrixXd;

double lx, ly;
int Nx, Ny;
int Nb, Ne, Nlb;
double h1, h2;
string shape, u0_str, f_str, delf_str;
string BCNodes, BCValues;
string u_AB, u_BC, u_DC, u_AD;
double abs_tol, rel_tol;
string output_file_dir;

void readFile(char* fileName) {
	xml_document<> doc;
	xml_node<> * root_node;
	// Read the xml file into a vector
    file<> theFile(fileName);
    doc.parse<0>(theFile.data());

	root_node = doc.first_node("MyConfig");
	// Iterate over the data
	for (xml_node<> * data_node = root_node->first_node("Data"); data_node; data_node = data_node->next_sibling()) {
        char* name = data_node->first_attribute("name")->value();
        char* value = data_node->first_attribute("value")->value();

        if (!strcmp(name, "lx")) lx = atof(value);
        else if (!strcmp(name, "ly")) ly = atof(value);
		else if (!strcmp(name, "Nx")) Nx = atoi(value);
		else if (!strcmp(name, "Ny")) Ny = atoi(value);
        else if (!strcmp(name, "u_AB")) u_AB = value;
        else if (!strcmp(name, "u_BC")) u_BC = value;
        else if (!strcmp(name, "u_DC")) u_DC = value;
        else if (!strcmp(name, "u_AD")) u_AD = value;
        else if (!strcmp(name, "BCNodes")) BCNodes = value;
        else if (!strcmp(name, "BCValues")) BCValues = value;
		else if (!strcmp(name, "shape")) shape = value;
		else if (!strcmp(name, "guess")) u0_str = value;
		else if (!strcmp(name, "source")) f_str = value;
		else if (!strcmp(name, "source_derivative")) delf_str = value;
		else if (!strcmp(name, "abs_tol")) abs_tol = atof(value);
		else if (!strcmp(name, "rel_tol")) rel_tol = atof(value);
		else if (!strcmp(name, "output_file_dir")) output_file_dir = value;
	}

    Nb = (Nx + 1) * (Ny + 1);
    Ne = Nx * Ny;
    Nlb = 4;
    h1 = lx / Nx;
    h2 = ly / Ny;
}

void setU_0(VectorXd& U_0) {
    double x,y;
    symbol_table symbols;
    symbols.add_variable("x", x);
    symbols.add_variable("y", y);

    expression u0;
    u0.register_symbol_table(symbols);

    parser Parser;
    Parser.compile(u0_str, u0);

    for (int n = 0; n < Nb; n++) {
        x = (n % (Nx + 1)) * h1;
        y = (n / (Nx + 1)) * h2;
        U_0(n) = u0.value();
    }
}

int main(int argc, char* argv[]) {
    readFile(argv[1]);

    if (shape == "rectangle") {
        VectorXd U_0(Nb);

        setU_0(U_0);

        MeshRect mesh = MeshRect(U_0);
        mesh.applyBCtoU();

        for (int i = 0; i < MAX_ITERATIONS; i++) {
            cout << "ITERATION " << i+1 << endl;

            mesh.calculateAF();
            mesh.applyBCtoDelU();

            mesh.addDelU();

            if (mesh.endConditionMet()){
                break;
            }
        }

        mesh.output();
    }
}