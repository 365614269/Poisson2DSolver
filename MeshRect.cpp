#include "include/MeshRect.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

MeshRect::MeshRect(VectorXd& U_) {
    this->Uv = U_;
    this->Fv = VectorXd::Zero(Nb);
    this->delUv = VectorXd::Zero(Nb);
    this->stiffness = MatrixXd::Zero(Nb, Nb);

    this->f_expr = this->parse_f();
    this->delf_expr = this->parse_delf();

    for (int i = 0; i < Ne; i++) {
        this->elements.push_back(Element(i));
    }

    for (int i = 0; i < Nb; i++) {
        int x = i / (Nx + 1);
        int y = i % (Nx + 1);

        if (y == 0) leftBNodes.push_back(i);
        else if (y == Nx) rightBNodes.push_back(i);
        if (x == 0) upperBNodes.push_back(i);
        else if (x == Ny) lowerBNodes.push_back(i);
    }
}

vector<string> MeshRect::split(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

Element& MeshRect::getElement(int n) {
    return this->elements[n];
}

int MeshRect::Tb(int localNodeIndex, int i, int j) {
    int offset[4][2] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
    int x = i + offset[localNodeIndex][0];
    int y = j + offset[localNodeIndex][1];

    return x * (Nx + 1) + y;  // exchange index
}

int MeshRect::Tb(int localNodeIndex, int n) {
    int offset[4][2] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
    pair<int, int> pos = make_pair(n / Nx, n % Nx);  // exchange index

    int x = pos.first + offset[localNodeIndex][0];
    int y = pos.second + offset[localNodeIndex][1];

    return x * (Nx + 1) + y;  // exchange index
}

void MeshRect::addAF(int elemIndex) {
    Element element = this->getElement(elemIndex);
    double ans1, ans2, ans3;
    double x1,x2,y1,y2;

    pair<int, int> index = make_pair(elemIndex / Nx, elemIndex % Nx);  // exchange index
    x1 = index.second * h1;
    x2 = (index.second + 1) * h1;
    y1 = (Ny - index.first - 1) * h2;
    y2 = (Ny - index.first) * h2;

    for (int alpha = 0; alpha < Nlb; alpha++) {
        ans3 = 0;

        Node& alphaNode = element.getNode(alpha);
        for (int beta = 0; beta < Nlb; beta++) {
            Node& betaNode = element.getNode(beta);

            ans1 = this->integrateDelpsi(betaNode, alphaNode, x1, x2, y1, y2);
            ans2 = this->integratePsiDelF(alphaNode, betaNode, x1, x2, y1, y2);
            ans3 += ans1 * this->Uv(Tb(beta, elemIndex));

            double ans = -ans1-ans2;
            this->stiffness(Tb(beta, elemIndex), Tb(alpha, elemIndex)) += ans;
        }

        ans3 += this->integrateFPsi(alphaNode, x1, x2, y1, y2);
        this->Fv(Tb(alpha, elemIndex)) += ans3;
    }
}

void MeshRect::calculateAF() {
    this->stiffness = MatrixXd::Zero(Nb, Nb);
    this->Fv = VectorXd::Zero(Nb);

    for (int i = 0; i < Ne; i++) {
        this->addAF(i);
    }
}

void MeshRect::applyBCtoU() {
    this->applyLeftBCtoU();
    this->applyRightBCtoU();
    this->applyUpperBCtoU();
    this->applyLowerBCtoU();
}

void MeshRect::applyUpperBCtoU() {
    if (u_AB == "") return;

    vector<string> cons = MeshRect::split(u_AB, ";");

    for (int i = 0; i < cons.size(); i++) {
        vector<string> con = MeshRect::split(cons[i], ",");

        int lb = stoi(con[0].substr(0, 3));
        int ub = stoi(con[0].substr(4, 3));
        double value = stof(con[1]);

        for (int j = 0; j < upperBNodes.size(); j++) {
            double percentage = ((j + 0.0) / upperBNodes.size()) * 100;
            int node = upperBNodes[j];

            if (percentage >= lb && percentage <= ub) {
                this->Uv(node) = value;
                this->BCNodes.push_back(node);
            }
        }
    }
}

void MeshRect::applyRightBCtoU() {
    if (u_BC == "") return;

    vector<string> cons = MeshRect::split(u_BC, ";");

    for (int i = 0; i < cons.size(); i++) {
        vector<string> con = MeshRect::split(cons[i], ",");
        int lb = stoi(con[0].substr(0, 3));
        int ub = stoi(con[0].substr(4, 3));
        double value = stof(con[1]);

        for (int j = 0; j < rightBNodes.size(); j++) {
            double percentage = ((j + 0.0) / rightBNodes.size()) * 100;
            int node = rightBNodes[j];

            if (percentage >= lb && percentage <= ub) {
                // cout << "RBNode " << j << " " << node << endl;
                this->Uv(node) = value;
                this->BCNodes.push_back(node);
            }
        }
    }
}

void MeshRect::applyLeftBCtoU() {
    if (u_AD == "") return;

    vector<string> cons = MeshRect::split(u_AD, ";");

    for (int i = 0; i < cons.size(); i++) {
        vector<string> con = MeshRect::split(cons[i], ",");
        int lb = stoi(con[0].substr(0, 3));
        int ub = stoi(con[0].substr(4, 3));
        double value = stof(con[1]);

        for (int j = 0; j < leftBNodes.size(); j++) {
            double percentage = ((j + 0.0) / leftBNodes.size()) * 100;
            int node = leftBNodes[j];

            if (percentage >= lb && percentage <= ub) {
                // cout << "LBNode " << node << endl;
                this->Uv(node) = value;
                this->BCNodes.push_back(node);
            }
        }
    }
}

void MeshRect::applyLowerBCtoU() {
    if (u_DC == "") return;

    vector<string> cons = MeshRect::split(u_DC, ";");

    for (int i = 0; i < cons.size(); i++) {
        vector<string> con = MeshRect::split(cons[i], ",");

        int lb = stoi(con[0].substr(0, 3));
        int ub = stoi(con[0].substr(4, 3));
        double value = stof(con[1]);

        for (int j = 0; j < lowerBNodes.size(); j++) {
            double percentage = ((j + 0.0) / lowerBNodes.size()) * 100;
            int node = lowerBNodes[j];

            if (percentage >= lb && percentage <= ub) {
                this->Uv(node) = value;
                this->BCNodes.push_back(node);
            }
        }
    }
}

void MeshRect::applyBCtoDelU() {
    for (int i = 0; i < BCNodes.size(); i++) {
        int node = BCNodes[i];

        // cout << node << endl;
        
        for (int i = 0; i < Nb; i++) {
            this->stiffness(node,i) = 0;
        }
        this->stiffness(node,node) = 1;
        this->Fv(node) = 0;
    }
}

void MeshRect::addDelU() {
    VectorXd delU = (this->stiffness).colPivHouseholderQr().solve(this->Fv);
    this->delUv = delU;
    this->Uv += delU;
}

bool MeshRect::endConditionMet() {
    double abs_err = this->delUv.norm() / this->Uv.norm();
    double rel_err = this->Fv.norm();
    cout << "Abs. Error: " << abs_err << " ;;; Rel. Error: " << rel_err << endl;

    return (abs_err < abs_tol) && (rel_err < rel_tol);
}

pair<int, int> MeshRect::elem(double x, double y) {
    int a = (ly - y) / h2;
    int b = x / h1;

    a = min(Ny, a);
    b = min(Nx, b);

    return make_pair(a, b);
}

double MeshRect::u(double x, double y){
    //find which element (x, y)is in.
    pair<int, int> index = this->elem(x, y);
    Element element = this->getElement(index.first * Nx + index.second);  // exchange index

    VectorXd v1(Nlb);
    VectorXd v2(Nlb);

    // add up the 4 components * shape functions of the 4 nodes

    for (int i = 0; i < Nlb; i++) {
        v1(i) = element.getNode(i).psi(x, y);
        int n = Tb(i, index.first, index.second);
        v2(i) = this->Uv(n);
    }

    return v1.dot(v2);
}

expression MeshRect::parse_f() {
    symbol_table symbols;
    symbols.add_variable("u", this->u_val);

    expression expr;
    expr.register_symbol_table(symbols);

    parser Parser;
    Parser.compile(f_str, expr);

    return expr;
}

expression MeshRect::parse_delf() {
    symbol_table symbols;
    symbols.add_variable("u", this->u_val);

    expression expr;
    expr.register_symbol_table(symbols);

    parser Parser;
    Parser.compile(delf_str, expr);

    return expr;
}

double MeshRect::f(double x, double y) {
    this->u_val = this->u(x, y);
    return f_expr.value();
}

double MeshRect::delf(double x, double y) {
    this->u_val = this->u(x, y);
    return delf_expr.value();
}

double MeshRect::integratePsiDelF(Node& node1, Node& node2, double a1, double b1, double a2, double b2) {
    double diff1 = (b1 - a1) / 2;
    double avg1 = (b1 + a1) / 2;
    double diff2 = (b2 - a2) / 2;
    double avg2 = (b2 + a2) / 2;

    double ans = 0;
    double x,y;

    for (int i = 0; i < PRECISION; i++) {
        for (int j = 0; j < PRECISION; j++) {
            x = diff1 * abscissae[i] + avg1;
            y = diff2 * abscissae[j] + avg2;
            ans += weights[i] * weights[j] * node1.psi(x, y) * node2.psi(x, y) * this->delf(x, y);
        }
    }

    return diff1 * diff2 * ans;
}

double MeshRect::integrateDelpsi(Node& node1, Node& node2, double a1, double b1, double a2, double b2) {
    double diff1 = (b1 - a1) / 2;
    double avg1 = (b1 + a1) / 2;
    double diff2 = (b2 - a2) / 2;
    double avg2 = (b2 + a2) / 2;

    double ans = 0;
    double x,y;

    for (int i = 0; i < PRECISION; i++) {
        for (int j = 0; j < PRECISION; j++) {
            x = diff1 * abscissae[i] + avg1;
            y = diff2 * abscissae[j] + avg2;
            ans += weights[i] * weights[j] * node1.delpsi(x, y).dot(node2.delpsi(x, y));
        }
    }

    return diff1 * diff2 * ans;
}

double MeshRect::integrateFPsi(Node& node, double a1, double b1, double a2, double b2) {
    double diff1 = (b1 - a1) / 2;
    double avg1 = (b1 + a1) / 2;
    double diff2 = (b2 - a2) / 2;
    double avg2 = (b2 + a2) / 2;

    double ans = 0;
    double x,y;

    for (int i = 0; i < PRECISION; i++) {
        for (int j = 0; j < PRECISION; j++) {
            x = diff1 * abscissae[i] + avg1;
            y = diff2 * abscissae[j] + avg2;            
            ans += weights[i] * weights[j] * node.psi(x, y) * this->f(x, y);
        }
    }

    return diff1 * diff2 * ans;
}

void MeshRect::output() {
    ofstream fout(output_file_dir);

    fout << "# vtk DataFile Version 3.0" << endl;
    fout << "2D Poisson Equation Numeric Solution" << endl;
    fout << "ASCII" << endl;
    fout << endl;
    fout << "DATASET UNSTRUCTURED_GRID" << endl;
    fout << "POINTS " << Nb << " float" << endl;

    for (int n = 0; n < Nb; n++) {
        fout << (n % (Nx + 1)) * h1 << ' ' << ly - (n / (Nx + 1)) * h2 << " 2" << endl;
    }

    fout << endl;
    fout << "CELLS " << Ne << " " << (Nlb + 1) * Ne << endl;

    for (int n = 0; n < Ne; n++) {
        fout << Nlb << " ";

        for (int i = 0; i < Nlb; i++) {
            fout << Tb(i, n) << " ";
        }

        fout << endl;
    }

    fout << endl;
    fout << "CELL_TYPES " << Ne << endl;
    
    for (int n = 0; n < Ne; n++) {
        fout << 9 << endl;
    }

    fout << endl;
    fout << "POINT_DATA " << Nb << endl;
    fout << "SCALARS temprerature float 1" << endl;
    fout << "LOOKUP_TABLE default" << endl;

    for (int n = 0; n < Nb; n++) {
        fout << this->Uv(n) << endl;
    }

    fout.close();
}