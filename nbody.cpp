#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
// #include "matplotlibcpp.h"


using namespace std;
// namespace plt = matplotlibcpp;
using doubles = std::vector<double>;

int clearMyPlot();
int writePlotToMyPlot(const doubles& xs);
doubles linspace(double start, double end, int numPoints);


const double G = 1;

/* MY TYPES */
struct Vector {
    double x, y, z;

    double dist() const {
        return sqrt(x*x + y*y + z*z);
    }
    Vector operator+(const Vector& other) const {
        return Vector{x+other.x, y+other.y, z+other.z};
    }
    Vector operator-(const Vector& other) const {
        return Vector{x-other.x, y-other.y, z-other.z};
    }
    Vector operator*(double scalar) const {
        return Vector{x*scalar, y*scalar, z*scalar};
    }
    friend Vector operator*(double scalar, const Vector& v) {
        return v*scalar;
    }
    Vector operator/(double scalar) const {
        return Vector{x/scalar, y/scalar, z/scalar};
    }
    Vector& operator+=(const Vector& other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }
    Vector& operator-=(const Vector& other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector& v) {
        os << "V (" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }
};
using Point = Vector;
using Vectors = std::vector<Vector>;


class Particle {
protected:
    static double totalMass;

    int id;
    double mass;
    Point p;
    Point v;
public: 
    Particle(int id, double mass, Point p, Point v) : id(id), mass(mass), p(p), v(v) {
        totalMass += mass;
    }

    friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
        os << "Particle " << p.id << " with mass " << p.mass << " at position (" << p.p.x << ", " << p.p.y << ", " << p.p.z << ") and velocity (" << p.v.x << ", " << p.v.y << ", " << p.v.z << ")";
        return os;
    }

    /* two Treeparticles would be equal with this definition! */
    bool operator==(const Particle& other) const {
        return id == other.id;
    }

    Vector getForceTo(const Particle& other, double eps) const {
        Vector r12 = other.p - p;
        double r = r12.dist();
        double force = G * mass * other.mass / (r*r + eps*eps);
        return force * r12/r;
    }

    double getDist() const {
        return p.dist();
    }
    double getMass() const {
        return mass;
    }
    Point getP() const {
        return p;
    }

    static double getTotalMass() {
        return totalMass;
    }

};
double Particle::totalMass = 0;
using Particles = std::vector<Particle>;

class Treeparticle : public Particle {
public: 
    Treeparticle() : Particle(-1, 0, Point{0, 0, 0}, Point{0, 0, 0}) {}
    Treeparticle(double mass, Vector p) : Particle(-1, mass, p, Point{0, 0, 0}) {}

    void add(const Treeparticle& other) {
        double newMass = mass + other.mass;
        if (newMass==0) return;
        Vector newP = (mass*p + other.mass*other.p) / newMass;
        mass = newMass;
        p = newP;
    }
};


Particles readIn(string filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    string line;
    getline(inputFile, line);           // skip first line
    // getline(inputFile, line);           // skip second line

    Particles particles;
    int i=0;
    while (getline(inputFile, line)) {
        istringstream iss(line);
        int id; iss >> id;
        double mass; iss >> mass;
        Point p; iss >> p.x >> p.y >> p.z;
        Point v; iss >> v.x >> v.y >> v.z;

        particles.push_back(Particle(id, mass, p, v));

        i++;
    }

    return particles;
}


/* TASK 1.1 */
doubles computeShells(const Particles& particles, double prec=0.1) {
    doubles shells;
    for (const auto& p : particles) {
        double dist = p.getDist();
        double shell = floor(dist/prec);
        if (shell >= shells.size()) {
            shells.resize(shell+1);
        }

        shells[shell] += p.getMass();
    }

    doubles densities(shells.size());
    for (uint i=0; i<shells.size(); i++) {
        double volume = 4./3. * M_PI * (pow((i+1)*prec, 3) - pow(i*prec, 3));
        densities[i] = shells[i] / volume;
    }
    return densities;
}
double rho(double r) {
    double M = Particle::getTotalMass();
    double a = .1;
    return M/(2*M_PI) * (a/r) * (1/(r+a)/(r+a)/(r+a));
}
vector<double> computeRhos(int size, double prec) {
    vector<double> rhos(size);
    for (int i=0; i<size; i++) {
        rhos[i] = rho(i*prec);
    }
    return rhos;
}
doubles computePoissonErrors(const doubles& shells, const doubles& lambdas) {
    /* in #of standard-deviations */
    doubles errors(shells.size());
    for (uint i=0; i<shells.size(); i++) {
        errors[i] = abs((shells[i] - lambdas[i]) / sqrt(lambdas[i]));
    }
    return errors;
}



/* TASK 1.2 */
double getHmNr(const Particles& particles) {
    /* calculate number of particles within the half-mass radius */
    if (particles[13].getMass() != particles[42].getMass()) throw std::runtime_error("Masses are not equal");
    
    double totalMass = Particle::getTotalMass();
    double centralMass = particles[0].getMass();
    double smallMass = particles[1].getMass();

    /* totalMass/2 = centralMass + i*smallMass */
    /* we need to return i+1 because we have center together with i small masses*/
    int halfMass_i = (totalMass/2. - centralMass) / smallMass + 1;
    return halfMass_i;
}
double getRhm(const Particles& particles) {
    /* calculate half-mass radius R_hm */
    int half = getHmNr(particles);
    // int half = particles.size()/2;          // can use half because all have equal mass
    doubles dists;
    for (const auto& p : particles) {
        dists.push_back(p.getDist());
    }
    nth_element(dists.begin(), dists.begin()+half-1, dists.end());
    return dists[half];
}

double getSoftening(const Particles& particles) {
    /* d = (V/N)^(1/3) */
    double Rhm = getRhm(particles);
    double V = 4/3 * M_PI * Rhm*Rhm*Rhm;
    int N = getHmNr(particles);

    return pow(V/N, 1/3);
}

Vectors getDirectForces(const Particles& particles, double softeningMultiplier=1) {
    double softening = getSoftening(particles)*softeningMultiplier;

    Vectors forces(particles.size());
    for (uint i=0; i<particles.size(); i++) {
        Particle p1 = particles[i];
        for (uint j=i+1; j<particles.size(); j++) {
            Particle p2 = particles[j];

            Vector force = p1.getForceTo(p2, softening);

            forces[i] += force;
            forces[j] -= force;
        }
    }
    return forces;
}

int getAnalyticalForce(const Particles& particles) {
    /* apply Newton II for spherical potentials */
    /* \vec F(r) = -GM(r)/r^3 * \vec r with M(r)=mass within radius r */
    return 0;
}

double getRelaxationTimescale(const Particles& particles) {
    /* t_relax = N/(8 ln N) * t_cross */
    double Rhm = getRhm(particles);
    double mass_halfMassRadius = Particle::getTotalMass() / 2;
    double v_c = sqrt(G*mass_halfMassRadius/Rhm);

    int N = particles.size();           // calculate with central obj.?
    double t_cross = 2*Rhm / v_c;       // calculate for half-mass radius ?

    double t_relax = N/(8*log(N)) * t_cross;

    return t_relax;
}

/* TASK 2 -- tree-code */
class Range {
private:
    double xmin, xmax;
    double ymin, ymax;
public:
    Range(double xmin, double xmax, double ymin, double ymax) : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax) {}
    bool contains(const Point& p) const {
        return xmin <= p.x && p.x <= xmax && ymin <= p.y && p.y <= ymax;
    }
    int getQuadrant(const Point& p) const {
        // could check if contains(p) first
        double xmid = (xmin + xmax) / 2;
        double ymid = (ymin + ymax) / 2;
        return (p.x > xmid) + 2*(p.y > ymid);
    }
    int getQuadrant(const Range& r) const {
        // could check if contains(r) first
        return (r.xmin != xmin) + 2*(r.ymin != ymin);
    }
    Range getChildRange(int quadrant) const {
        double xmid = (xmin + xmax) / 2;
        double ymid = (ymin + ymax) / 2;
        switch (quadrant) {
            case 0: return Range(xmin, xmid, ymin, ymid);
            case 1: return Range(xmid, xmax, ymin, ymid);
            case 2: return Range(xmin, xmid, ymid, ymax);
            case 3: return Range(xmid, xmax, ymid, ymax);
            default: throw std::runtime_error("Invalid quadrant");
        }
    }
    double openingAngle(Point p) const {
        if (contains(p)) return 2*M_PI;
        double poss[6][2][2] = {
            {{xmin, ymin}, {xmin, ymax}},
            {{xmin, ymax}, {xmax, ymax}},
            {{xmax, ymax}, {xmax, ymin}},
            {{xmax, ymin}, {xmin, ymin}},
            {{xmin, ymin}, {xmax, ymax}},
            {{xmax, ymin}, {xmin, ymax}}
        };
        double max = 0;
        for (const auto& pos : poss) {
            double x1 = pos[0][0], y1 = pos[0][1];
            double x2 = pos[1][0], y2 = pos[1][1];

            double arc1 = atan2(p.y-y1, p.x-x1);
            double arc2 = atan2(p.y-y2, p.x-x2);
            double arc = min(abs(arc1-arc2), 2*M_PI - abs(arc1-arc2));

            if (arc>max) max = arc;
        }
        return max;
    }

    friend std::ostream& operator<<(std::ostream& os, const Range& r) {
        os << "[(" << r.xmin << "," << r.xmax << ") x (" << r.ymin << "," << r.ymax << ")]";
        return os;
    }
};

/* 
|-----|-----|
|  2  |  3  |
|-----|-----|
|  0  |  1  |
|-----|-----|
*/
enum class NodeType { LEAF, NODE };
class Node;
class Leaf; 

class Tree {
protected:
    Range range;
    Node* parent;
    // Point centerOfMass;
    // double mass;
public:
    Treeparticle treeParticle;
    Tree(const Range& range, Node* parent) : range(range), parent(parent) {}

    virtual ~Tree() = default;
    virtual bool add(const Particle& p) = 0;
    virtual NodeType getType() const = 0;
    virtual void print(std::ostream& os) const = 0;

    friend std::ostream& operator<<(std::ostream&os, const Tree& tree) {
        tree.print(os);
        return os;
    }

    Range getRange() const {
        return range;
    }

    // void setTreeAsParticle(double mass, Point centerOfMass) {
    //     treeParticle = Treeparticle(mass, centerOfMass);
    //     // this->mass = mass;
    //     // this->centerOfMass = centerOfMass;
    // }

    // Treeparticle getTreeAsParticle() const {
    //     // return Particle(-1, mass, centerOfMass, Point{0, 0, 0});
    //     return treeParticle;
    // }

    static Tree* makeRoot(const Range& range);
    // static Tree* makeRoot(const Range& range) {
    //     return new Node(range, nullptr);
    // }
};

class Leaf : public Tree {
private:
    vector<const Particle*> ps;
public:
    Leaf(const Range& range, Node* parent) : Tree(range, parent) {}

    bool add(const Particle&);
    // bool add(const Particle& p) {
    //     ps.push_back(&p);
    //     if (ps.size() > 8) {
    //         /* upgrade from Leaf to Node */
    //         Node* newNode = new Node(range, parent);
    //         parent->setNewChild(newNode);
    //         for (const auto& p : ps) {
    //             newNode->add(*p);
    //         }

    //         delete this;                    // dangerous
    //     }
    //     return true;
    // }

    NodeType getType() const override {
        return NodeType::LEAF;
    }

    void print(std::ostream& os) const override {
        os << "Leaf with " << getNrParticles() << " particles in range " << range;
    }

    vector<const Particle*> getParticles() const {
        return ps;
    }

    int getNrParticles() const {
        return ps.size();
    }
};

class Node : public Tree {
protected:
    Tree* children[4];
public:
    Node(const Range& range, Node* parent) : Tree(range, parent) {
        for (uint i=0; i<4; i++) {
            children[i] = new Leaf(range.getChildRange(i), this);
        }
    }
    ~Node() override {
        for (auto& child : children) {
            delete child;
        }
    }

    bool add(const Particle& p) {
        int q = range.getQuadrant(p.getP());
        return children[q]->add(p);
    }

    NodeType getType() const override {
        return NodeType::NODE;
    }

    void print(std::ostream& os) const override {
        os << "Node with range " << range;
    }

    bool setNewChild(Node* child) {
        int q = range.getQuadrant(child->range);
        children[q] = child;
        return true;
    }

    Tree* getChild(int quadrant) const {
        return children[quadrant];
    }
};
Tree* Tree::makeRoot(const Range& range) {
        return new Node(range, nullptr);
}

bool Leaf::add(const Particle& p) {
    ps.push_back(&p);
    if (ps.size() > 8) {
        /* upgrade from Leaf to Node */
        Node* newNode = new Node(range, parent);
        parent->setNewChild(newNode);
        for (const auto& p : ps) {
            newNode->add(*p);
        }

        delete this;                    // dangerous
    }
    return true;
}

/* for seeing if the tree was built correctly */
int treeWalk(const Tree* tree) {
    if (tree->getType() == NodeType::LEAF) {
        return dynamic_cast<const Leaf*>(tree)->getNrParticles();
    } else {
        int sum = 0;
        for (int i=0; i<4; i++) {
            sum += treeWalk(dynamic_cast<const Node*>(tree)->getChild(i));
        }
        if (sum > 220) {
            std::cout << "Treenode with " << sum << " particles: " << *tree << endl;
        }
        return sum;
    }
}
Tree* getTree(const Particles& particles) {
    Range range(-1, 1, -1, 1);
    for (const auto& p : particles) {
        if (!range.contains(p.getP())) throw std::runtime_error("Particle not in range");
    }

    cout << "Creating tree..." << endl;
    int i=0;

    Tree* root = Tree::makeRoot(range);
    for (const auto& p : particles) {
        root->add(p);
        i++;
    }
    cout << "Added all particles..." << endl;
    return root;
}

Treeparticle computeTreeCodeMasses(Tree* tree) {                // could be made const again and put into tree class??
    if (tree->getType() == NodeType::LEAF) {
        Leaf* leaf = dynamic_cast<Leaf*>(tree);
        for (const auto& p : leaf->getParticles()) {
            Treeparticle tp = Treeparticle(p->getMass(), p->getP());
            leaf->treeParticle.add(tp);
        }
    } 
    else {
        Node* node = dynamic_cast<Node*>(tree);
        for (int i=0; i<4; i++) {
            Treeparticle tp = computeTreeCodeMasses(node->getChild(i));
            node->treeParticle.add(tp);
        }
    }
    return tree->treeParticle;
}

Vector getForceRec(const Particle& p, const Tree* tree, double angle, double eps) {
    if (tree->getRange().openingAngle(p.getP()) > angle) {
        Vector force{0,0,0};
        if (tree->getType() == NodeType::LEAF) {
            const Leaf* leaf = dynamic_cast<const Leaf*>(tree);
            for (const auto& p2 : leaf->getParticles()) {
                if (p == *p2) continue;
                force += p.getForceTo(*p2, eps);
            }
        } else {
            const Node* node = dynamic_cast<const Node*>(tree);
            for (int i=0; i<4; i++) {
                force += getForceRec(p, node->getChild(i), angle, eps);
                // std::cout << "Force from particle " << p << " to range " << node->getChild(i)->getRange() << ": " << force << endl;
            }
        }
        return force;
    }
    else {
        return p.getForceTo(tree->treeParticle, eps);
    }

}

Vectors getTreeCodeForces(const Particles& particles, const Tree* tree, double softeningMultiplier=1, double angleMultiplier=1) {
    double softening = getSoftening(particles)*softeningMultiplier;
    double angle = M_PI/4 * angleMultiplier; // 45 degrees
    Vectors forces;

    for (const auto& p : particles) {
        forces.push_back(getForceRec(p, tree, angle, softening));
    }
        
    return forces;
}

void compareForces(const Vectors& v1, const Vectors& v2) {
    if (v1.size() != v2.size()) throw std::runtime_error("Vectors have different sizes");
    double sum = 0;
    for (uint i=0; i<v1.size(); i++) {
        Vector diff = v1[i] - v2[i];
        // if (diff.dist() > 1e-2) {
        //     std::cout << "Difference in force for particle " << i << ": " << diff << endl;
        // }
        sum += diff.dist();
    }
    std::cout << "Total difference in forces: " << sum << endl;
}






int main() {
    /*
    - central mass is 1 in both datasets
    - total particle mass is 10 in both datasets
    */

    string dataset = "data0.txt";
    Particles particles = readIn(dataset);

    /* TESTING DATA INTEGRITY */
    std::cout << "Total mass: " << Particle::getTotalMass() << endl;
    std::cout << particles[0] << endl;
    std::cout << "--------------------" << endl;

    /* TASK 1.1 */
    // double prec = .06;
    // doubles shells = computeShells(particles, prec);
    // doubles rhos = computeRhos(shells.size(), prec);
    // doubles xs = linspace(0, prec*shells.size(), shells.size());
    // doubles poissonErrors = computePoissonErrors(shells, rhos);
    
    // clearMyPlot();
    // writePlotToMyPlot(xs);
    // writePlotToMyPlot(shells);
    // writePlotToMyPlot(rhos);
    // writePlotToMyPlot(poissonErrors);

    /* TASK 1.2 */
    Vectors v1 = getDirectForces(particles, 0);
    // Vectors v2 = getDirectForces(particles, .1);
    Vectors v3 = getDirectForces(particles, 1);
    // Vectors v4 = getDirectForces(particles, 10);

    /* TODO: IMPLEMENT ANALYTICAL FORCE STUFF */

    // double relaxationTime = getRelaxationTimescale(particles);
    // std::cout << "Relaxation time: " << relaxationTime << endl;
    // std::cout << "For an increased softening, one can expect a longer relaxation time,\n"
    // "as the particles are less affected by each other (and the relaxation time gives an\n"
    // "estimate of how long it takes for the particles to change their velocities due to\n"
    // "gravitational interactions with other particles)." << endl;

    /* TASK 2 -- tree-code */
    Tree* root = getTree(particles);        // needs to be deleted!
    // treeWalk(root);
    computeTreeCodeMasses(root);
    // std::cout << "Tree particle of root: " << root->treeParticle << endl;

    Vectors v5 = getTreeCodeForces(particles, root, 0, 1);
    Vectors v5_ = getTreeCodeForces(particles, root, 0, .3);
    Vectors v5__ = getTreeCodeForces(particles, root, 0, .1);
    Vectors v6 = getTreeCodeForces(particles, root, .1, 1);
    Vectors v7 = getTreeCodeForces(particles, root, 1, 1);
    Vectors v7_ = getTreeCodeForces(particles, root, 1, .1);
    Vectors v8 = getTreeCodeForces(particles, root, 10, 1);

    std::cout << "Comparing forces for eps=0,ang= 1 "; compareForces(v1, v5);
    std::cout << "Comparing forces for eps=0,ang=.3 "; compareForces(v1, v5_);
    std::cout << "Comparing forces for eps=0,ang=.1 "; compareForces(v1, v5__);
    std::cout << "Comparing forces for eps=1,ang= 1 "; compareForces(v3, v7);
    std::cout << "Comparing forces for eps=1,ang=.1 "; compareForces(v3, v7_);




    delete root;
    return 0;
}

int clearMyPlot() {
    std::ofstream clearFile("myplot.txt");
    if (!clearFile.is_open()) throw std::runtime_error("Could not open file: myplot.txt");
    clearFile << "";

    return 0;
}

int writePlotToMyPlot(const doubles& xs) {
    std::ofstream appendFile("myplot.txt", std::ios_base::app);
    if (!appendFile.is_open()) throw std::runtime_error("Could not open file: myplot.txt");
    for (double x : xs) appendFile << x << " ";
    appendFile << endl;

    return 0;
}

vector<double> linspace(double start, double end, int numPoints) {
    vector<double> xs(numPoints);
    double step = (end-start)/(numPoints-1);
    for (int i=0; i<numPoints; i++) xs[i] = start + i*step;
    return xs;
}