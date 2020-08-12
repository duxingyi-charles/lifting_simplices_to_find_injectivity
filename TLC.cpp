// 2D/3D lifted formulation
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <chrono>

#include <nlopt.hpp>



bool importData(const char* filename,
	std::vector<std::vector<double> > &restV,
	std::vector<std::vector<double> > &initV,
	std::vector<std::vector<unsigned> > &F,
	std::vector<unsigned> &handles
	)
{
	std::ifstream in_file(filename);

	if (! in_file.is_open()) {
        std::cerr << "Failed to open " << filename << "!" << std::endl;
        return false;
    }

    //read the file
    size_t n, ndim;
    // restV
    in_file >> n >> ndim;
    restV.resize(n);
    for (int i = 0; i < n; ++i)
    {
    	std::vector<double> v(ndim);
    	for (int j = 0; j < ndim; ++j)
    	{
    		in_file >> v[j];
    	}
    	restV[i] = v;
    }

    //initV
    in_file >> n >> ndim;
    initV.resize(n);
    for (int i = 0; i < n; ++i)
    {	
    	std::vector<double> v(ndim);
    	for (int j = 0; j < ndim; ++j)
    	{
    		in_file >> v[j];
    	}
    	initV[i] = v;
    }

    //F
    size_t simplexSize;
    in_file >> n >> simplexSize;
    F.resize(n);
    for (int i = 0; i < n; ++i)
    {	
    	std::vector<unsigned> v(simplexSize);
    	for (int j = 0; j < simplexSize; ++j)
    	{
    		in_file >> v[j];
    	}
    	F[i] = v;
    }

    //handles
    in_file >> n;
    handles.resize(n);
    for (int i = 0; i < n; ++i)
    {		
    	unsigned v;
    	in_file >> v;
    	handles[i] = v;
    }

    in_file.close();

    return true;
}

// solver options
class NloptOptionManager
{
public:
	//default options
	NloptOptionManager():
	form("Tutte"), alphaRatio(1e-6), alpha(-1),
	ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8),
	maxeval(10000), algorithm("LBFGS"), stopCode("all_good"), record()
	{};
	//import options from file
	explicit NloptOptionManager(const char* filename):
	form("Tutte"), alphaRatio(1e-6), alpha(-1),
	ftol_abs(1e-8), ftol_rel(1e-8), xtol_abs(1e-8), xtol_rel(1e-8),
	maxeval(10000), algorithm("LBFGS"), stopCode("all_good"), record()
	{
		if (!importOptions(filename))
		{
			std::cout << "NloptOptionManager Warn: default options are used." << std::endl;

		}
	};

	~NloptOptionManager() = default;

	// energy formulation options
	std::string form;
	double alphaRatio;
	double alpha;

	// optimization options
	double ftol_abs;
	double ftol_rel;
	double xtol_abs;
	double xtol_rel;
	int maxeval;
	std::string algorithm;
	std::string stopCode;
	std::vector<std::string> record;

	void printOptions()
	{
        std::cout << "form:\t" << form << "\n";
        std::cout << "alphaRatio:\t" << alphaRatio << "\n";
        std::cout << "alpha:\t"    <<  alpha    << "\n";
		std::cout << "ftol_abs:\t" <<  ftol_abs << "\n";
		std::cout << "ftol_rel:\t" <<  ftol_rel << "\n";
		std::cout << "xtol_abs:\t" <<  xtol_abs << "\n";
		std::cout << "xtol_rel:\t" <<  xtol_rel << "\n";
		std::cout << "maxeval:\t" <<  maxeval << "\n";
		std::cout << "algorithm:\t" <<  algorithm << "\n";
		std::cout << "stopCode:\t" <<  stopCode << "\n";
		std::cout << "record:  \t" <<  "{ ";
		for (std::vector<std::string>::iterator i = record.begin(); i != record.end(); ++i)
		{
			std::cout << *i << " ";
		}
		std::cout << "}" << std::endl;
	}

	bool importOptions(const char* filename)
	{
		//open the data file
		std::ifstream in_file(filename);

		if (! in_file.is_open()) 
		{
			std::cerr << "Failed to open solver option file:" << filename << "!" << std::endl;
			return false;
		}

		//read the file
		unsigned abnormal = 0;
		std::string optName;
		while(true)
		{
		    in_file >> optName;
            if (optName != "form")
            {
                abnormal = 1;
                break;
            }
            in_file >> form;

            in_file >> optName;
            if (optName != "alphaRatio")
            {
                abnormal = 2;
                break;
            }
            in_file >> alphaRatio;

            in_file >> optName;
            if (optName != "alpha")
            {
                abnormal = 3;
                break;
            }
            in_file >> alpha;

			in_file >> optName;
			if (optName != "ftol_abs")
			{
				abnormal = 4;
				break;
			}
			in_file >> ftol_abs;

			in_file >> optName;
			if (optName != "ftol_rel")
			{
				abnormal = 5;
				break;
			}
			in_file >> ftol_rel;

			in_file >> optName;
			if (optName != "xtol_abs")
			{
				abnormal = 6;
				break;
			}
			in_file >> xtol_abs;

			in_file >> optName;
			if (optName != "xtol_rel")
			{
				abnormal = 7;
				break;
			}
			in_file >> xtol_rel;

			in_file >> optName;
			if (optName != "algorithm")
			{
				abnormal = 8;
				break;
			}
			in_file >> algorithm;

			in_file >> optName;
			if (optName != "maxeval")
			{
				abnormal = 9;
				break;
			}
			in_file >> maxeval;

			in_file >> optName;
			if (optName != "stopCode")
			{
				abnormal = 10;
				break;
			}
			in_file >> stopCode;

			// record values
			in_file >> optName;
			if (optName != "record")
			{
				abnormal = 11;
				break;
			}

            in_file >> optName;
            if (optName != "vert") {
                abnormal = 12;
                break;
            }
            int selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("vert");
            }

            in_file >> optName;
            if (optName != "energy") {
                abnormal = 13;
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("energy");
            }

            in_file >> optName;
            if (optName != "minArea") {
                abnormal = 14;
                break;
            }
            selected = 0;
            in_file >> selected;
            if (selected > 0) {
                record.emplace_back("minArea");
            }

			break;
		}

		in_file.close();

		if (abnormal!=0)
		{
			std::cout << "Err:(" << abnormal << ") fail to import options from file. Check file format." << std::endl;
			return false;
		}

		return true;
	}
	
};


//helper

double distance2(const std::vector<double> &p1, const std::vector<double> &p2)
{
	unsigned n = p1.size();
	double dist = 0;
	for (int i = 0; i < n; ++i)
	{
		dist += (p1[i]-p2[i])*(p1[i]-p2[i]);
	}
	return dist;
}

// v1 - v2
std::vector<double> vec_substract(const std::vector<double> &v1, const std::vector<double> &v2)
{
	std::vector<double> v(v1.size());
	for (int i = 0; i < v1.size(); ++i)
	{
		v[i] = v1[i] - v2[i];
	}
	return v;
}

// v1 = v1 + v2
void vec_add_to(std::vector<double> &v1, const std::vector<double> &v2)
{
	for (int i = 0; i < v1.size(); ++i)
	{
		v1[i] += v2[i];
	}
}

// v = s*v
void vec_scale(std::vector<double> &v, double s)
{
	for (int i = 0; i < v.size(); ++i)
	{
		v[i] *= s;
	}
}

// v1 = v1 + s*v2
void vec_scale_add_to(std::vector<double> &v1, const std::vector<double> &v2, double s)
{
	for (int i = 0; i < v1.size(); ++i)
	{
		v1[i] += (s*v2[i]);
	}
}


// add m2 to m1
void mat_add_to(std::vector<std::vector<double> > &m1, const std::vector<std::vector<double> > &m2)
{
	unsigned m = m1.size();
	unsigned n = m1[1].size();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			m1[i][j] += m2[i][j];
		}
	}
}

// pts is a vector of 2D coordinates of 3 points of triangle
// return: 2 times of signed area of the triangle
double tri_signed_area(const std::vector<std::vector<double> > &pts)
{
	double x1=pts[0][0], y1=pts[0][1];
	double x2=pts[1][0], y2=pts[1][1];
	double x3=pts[2][0], y3=pts[2][1];
	double area = x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3);
	return area;
}

// pts is a vector of 3D coordinates of 4 points of a tet
// return: 6 times of signd volume of the tet
double tet_signed_volume(const std::vector<std::vector<double> > &pts)
{
	double x1=pts[0][0], y1=pts[0][1], z1=pts[0][2];
	double x2=pts[1][0], y2=pts[1][1], z2=pts[1][2];
	double x3=pts[2][0], y3=pts[2][1], z3=pts[2][2];
	double x4=pts[3][0], y4=pts[3][1], z4=pts[3][2];
	double vol = x4*(y3*(-z1 + z2) + y2*(z1 - z3) + y1*(-z2 + z3)) + 
   	x3*(y4*(z1 - z2) + y1*(z2 - z4) + y2*(-z1 + z4)) + 
   	x1*(y4*(z2 - z3) + y2*(z3 - z4) + y3*(-z2 + z4)) + 
   	x2*(y4*(-z1 + z3) + y3*(z1 - z4) + y1*(-z3 + z4));
   	return vol;
}

// input: 2D coordinates of 3 points of a triangle
// return: 2 times of signed volume of the triangle
double tri_signed_area(const std::vector<double> &p1,
	const std::vector<double> &p2,
	const std::vector<double> &p3)
{
	double x1=p1[0], y1=p1[1];
	double x2=p2[0], y2=p2[1];
	double x3=p3[0], y3=p3[1];
	double area = x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3);
	return area;
}

// input: 3D coordinates of 4 points of a tet
// return: 6 times of signed volume of the tet
double tet_signed_volume(const std::vector<double> &p1,
	const std::vector<double> &p2,
	const std::vector<double> &p3,
	const std::vector<double> &p4)
{
	double x1=p1[0], y1=p1[1], z1=p1[2];
	double x2=p2[0], y2=p2[1], z2=p2[2];
	double x3=p3[0], y3=p3[1], z3=p3[2];
	double x4=p4[0], y4=p4[1], z4=p4[2];
	double vol = x4*(y3*(-z1 + z2) + y2*(z1 - z3) + y1*(-z2 + z3)) + 
   	x3*(y4*(z1 - z2) + y1*(z2 - z4) + y2*(-z1 + z4)) + 
   	x1*(y4*(z2 - z3) + y2*(z3 - z4) + y3*(-z2 + z4)) + 
   	x2*(y4*(-z1 + z3) + y3*(z1 - z4) + y1*(-z3 + z4));
   	return vol;
}

// input: tri-mesh with vertices V and faces F
// return: total signed area of the tri-mesh
double total_tri_signed_area(const std::vector<std::vector<double> > &V,
                             const std::vector<std::vector<unsigned> > &F)
{
    unsigned nF = F.size();

    double area = 0;
    for (int i = 0; i < nF; ++i) {
        area += tri_signed_area(V[F[i][0]],V[F[i][1]],V[F[i][2]]);
    }
    area /= 2;

    return area;
}

// input: tet-mesh with vertices V and faces F
// return: total signed volume of the tet-mesh
double total_tet_signed_volume(const std::vector<std::vector<double> > &V,
                             const std::vector<std::vector<unsigned> > &F)
{
    unsigned nF = F.size();

    double area = 0;
    for (int i = 0; i < nF; ++i) {
        area += tet_signed_volume(V[F[i][0]],V[F[i][1]],V[F[i][2]],V[F[i][3]]);
    }
    area /= 6;

    return area;
}

// input: d is a vector of squared length of 3 edges of a triangle
// return: 4 times of the triangle area
double tri_area(const std::vector<double> &d)
{
	// more numerical robust Heron's formula
	// sort d1,d2,d3 as a >= b >= c
	double d1=d[0],d2=d[1],d3=d[2];
	double a,b,c;
	if (d1 > d2) { a = d1; b = d2; }
	else { a = d2; b = d1; }
	c = d3;
	if (d3 > b) {
	   c = b;
	   b = d3;
	   if (d3 > a) {
	       b = a;
	       a = d3;
	   }
	}

	a = sqrt(a);
	b = sqrt(b);
	c = sqrt(c);

    return sqrt(abs((a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))));
}

// input: d is a vector of squared length of 6 edges of a tet
// return: 12 times of the tet volume
double tet_volume(const std::vector<double> &d)
{
	double d0=d[0],d1=d[1],d2=d[2],d3=d[3],d4=d[4],d5=d[5];
	double det = -(d1*d1*d4) - d0*d0*d5 - d3*(d2*d2 + d2*(d3 - d4 - d5) + d4*d5) + 
	d1*(d2*(d3 + d4 - d5) + d4*(d3 - d4 + d5)) + 
	d0*((d3 + d4 - d5)*d5 + d2*(d3 - d4 + d5) + d1*(-d3 + d4 + d5));
	return sqrt(det);
}

// input: mesh with vertices V and cells F
// return (reference): vector of squared edge lengths
void compute_mesh_squared_edge_length(const std::vector<std::vector<double> > &V,
        const std::vector<std::vector<unsigned> > &F,
        std::vector<std::vector<double> > &edgeLen2)
{
    unsigned nF = F.size();
    unsigned simplexSize = F[0].size();
    unsigned nEdge = simplexSize * (simplexSize -1) /2;

    // compute squared edge length
    edgeLen2.resize(0);
    for (int i = 0; i < nF; ++i) {
        std::vector<double> Di;
        Di.reserve(nEdge);
        for (int j = 0; j < simplexSize; ++j)
        {
            for (int k = j+1; k < simplexSize; ++k)
            {
                Di.push_back(distance2(V[F[i][j]],V[F[i][k]]));
            }
        }
        edgeLen2.push_back(Di);
    }
}

// input: a tri-mesh with vertices V and faces F
// output: total unsigned area of the mesh
double total_tri_area(const std::vector<std::vector<double> > &V,
        const std::vector<std::vector<unsigned> > &F)
{
    unsigned nF = F.size();

    std::vector<std::vector<double> > edgeLen2;
    compute_mesh_squared_edge_length(V,F,edgeLen2);

    // compute triangle areas
    double area = 0;
    for (int i = 0; i < nF; ++i) {
        area += tri_area(edgeLen2[i]);
    }
    area /= 4;

    return  area;
}

// input: a tet-mesh with vertices V and tets F
// output: total unsigned volume of the mesh
double total_tet_volume(const std::vector<std::vector<double> > &V,
                      const std::vector<std::vector<unsigned> > &F)
{
    unsigned nF = F.size();

    std::vector<std::vector<double> > edgeLen2;
    compute_mesh_squared_edge_length(V,F,edgeLen2);

    // compute tet volumes
    double area = 0;
    for (int i = 0; i < nF; ++i) {
        area += tet_volume(edgeLen2[i]);
    }
    area /= 12;

    return  area;
}

std::vector<double> tri_grad(const std::vector<double> &d)
{
	double d0=d[0],d1=d[1],d2=d[2];
	std::vector<double> g(3);
	g[0] = 2*(d1+d2-d0);
	g[1] = 2*(d0+d2-d1);
	g[2] = 2*(d0+d1-d2);
	return g;
}

std::vector<double> tet_grad(const std::vector<double> &d)
{
	double d0=d[0],d1=d[1],d2=d[2],d3=d[3],d4=d[4],d5=d[5];
	std::vector<double> g(6);
	g[0] = (-d1 + d2)*(d3 - d4) + (-2.*d0 + d1 + d2 + d3 + d4 - d5)*d5;
	g[1] = (-d0 + d2)*(d3 - d5) + d4*(d0 - 2.*d1 + d2 + d3 - d4 + d5);
	g[2] = (-d0 + d1)*(d4 - d5) + d3*(d0 + d1 - 2.*d2 - d3 + d4 + d5);
	g[3] = (-d0 + d4)*(d1 - d5) + d2*(d0 + d1 - d2 - 2.*d3 + d4 + d5);
	g[4] = (-d0 + d3)*(d2 - d5) + d1*(d0 - d1 + d2 + d3 - 2.*d4 + d5);
	g[5] = (-d1 + d3)*(d2 - d4) + d0*(-d0 + d1 + d2 + d3 + d4 - 2.*d5);
	return g;
}

// 
class LiftedData
{
public:
	 LiftedData(std::vector<std::vector<double> > &restV,
	 	std::vector<std::vector<double> > &initV,
	 	std::vector<std::vector<unsigned> > &restF,
	 	std::vector<unsigned> &handles,
	 	std::string form,
	 	double alphaRatio,
	 	double alpha) :
	 V(initV), F(restF), freeI(0), solutionFound(false),
	 lastFunctionValue(HUGE_VAL), stopCode("none"),
	 nb_feval(0),nb_geval(0),
	 record_vert(false), record_energy(false), record_minArea(false),
	 vertRecord(0),energyRecord(0), minAreaRecord(0)
	 {
	 	// compute free indices
	 	unsigned nV = V.size();
	 	std::vector<bool> freeQ(nV, true);
	 	for (std::vector<unsigned>::iterator i = handles.begin(); i != handles.end(); ++i)
	 	{
	 		freeQ[*i] = false;
	 	}
	 	for (int i = 0; i < nV; ++i)
	 	{
	 		if (freeQ[i])
	 		{
	 			freeI.push_back(i);
	 		}
	 	}

	 	// compute alpha if it's not explicitly given
	 	double alphaFinal;
	 	if (alpha >=0) alphaFinal = alpha;
	 	else {  //alpha < 0. Deduce alpha from alphaRatio
	 	    alphaFinal = computeAlpha(restV,initV,F,form,alphaRatio);
	 	}
	 	std::cout << "alphaRatio: " << alphaRatio << std::endl;
	 	std::cout << "alpha: " << alphaFinal << std::endl;

	 	//compute squared edge length of rest mesh
	 	unsigned nF = F.size();
	 	unsigned simplexSize = F[0].size();
	 	unsigned nEdge;
	 	double a;
	 	if (simplexSize == 3) { //triangle
	 		nEdge = 3;
	 		a = alphaFinal;
	 	} 
	 	else { //tet
	 		nEdge = 6;
	 		a = cbrt(alphaFinal);
	 		a = a*a;
	 	}
	 	
	 	if (form == "harmonic")
	 	{
	 		for (int i = 0; i < nF; ++i)
	 		{
	 			std::vector<double> Di;
	 			Di.reserve(nEdge);
	 			for (int j = 0; j < simplexSize; ++j)
	 			{
	 				for (int k = j+1; k < simplexSize; ++k)
	 				{
	 					Di.push_back(a * distance2(restV[F[i][j]],restV[F[i][k]]));
	 				}
	 			}
	 			restD.push_back(Di);
	 		}
	 	}
	 	else  //Tutte
	 	{
	 		for (int i = 0; i < nF; ++i)
	 		{
	 			std::vector<double> Di(nEdge,a);
	 			restD.push_back(Di);
	 		}
	 	}
		

		//compute x0 from initV
		unsigned ndim = V[0].size();
		x0.resize(ndim*freeI.size());
		for (int i = 0; i < freeI.size(); ++i)
		{
			for (int j = 0; j < ndim; ++j)
			{
				x0[ndim*i+j] = V[freeI[i]][j];
			}
		}

	 };

	~ LiftedData() = default;
	std::vector<std::vector<double> > V;
	std::vector<unsigned> freeI;
	std::vector<std::vector<unsigned> > F;
	std::vector<std::vector<double> > restD;
	std::vector<double> x0;
	bool solutionFound;
	double lastFunctionValue;

	int nb_feval;
	int nb_geval;

	//from options
	std::string stopCode;
	//record data
	bool record_vert;
	bool record_energy;
	bool record_minArea;
	std::vector<std::vector<std::vector<double> > > vertRecord;
	std::vector<double> minAreaRecord;
	std::vector<double> energyRecord;

	// compute alpha from alphaRatio
    static double computeAlpha(const std::vector<std::vector<double> > &restV,
            const std::vector<std::vector<double> > &initV,
            const std::vector<std::vector<unsigned> > &Faces,
            std::string form, double alphaRatio)
    {
        unsigned nF = Faces.size();
        unsigned simplex_size = Faces[0].size();
        double rest_measure = 0;
        double init_measure = 0;
        if (simplex_size == 3) { // tri
            init_measure = total_tri_signed_area(initV,Faces);
            if (form == "harmonic") {
                rest_measure = total_tri_area(restV,Faces);
            }
            else { // Tutte form
                rest_measure = nF * sqrt(3) / 4;
            }
        }
        else { // tet
            init_measure = total_tet_signed_volume(initV,Faces);
            if (form == "harmonic") {
                rest_measure = total_tet_volume(restV,Faces);
            }
            else { // Tutte form
                rest_measure = nF * sqrt(2) / 12;
            }
        }

        // alpha
        return alphaRatio * init_measure / rest_measure;
    }


	// record information we cared about
	void set_record_flags(const std::vector<std::string> &record)
	{
		for (std::vector<std::string>::const_iterator i = record.begin(); i != record.end(); ++i)
		{
			if (*i == "vert")    record_vert = true;
			if (*i == "energy")  record_energy = true;
			if (*i == "minArea") record_minArea = true;
		}
	}

	void record()
	{
		if (record_vert) vertRecord.push_back(V);
		//you need to make sure lastFunctionValue is up-to-date
		if (record_energy) energyRecord.push_back(lastFunctionValue);
		if (record_minArea)
		{
			double minA = HUGE_VAL;
			double curA;
			unsigned nF = F.size();
			unsigned simplexSize = F[0].size();
			std::vector<std::vector<double> > pts(simplexSize);
			if (simplexSize == 3)
			{ //triMesh
				for (int i = 0; i < nF; ++i)
				{
					for (int j = 0; j < simplexSize; ++j)
					{
						pts[j] = V[F[i][j]];
					}
					curA = tri_signed_area(pts);
					minA = (curA < minA) ? curA : minA;
				}
				//tri_signed_area returns 2 times of the area
				minA = minA/2;
			}
			else { //tetMesh 
				for (int i = 0; i < nF; ++i)
				{
					for (int j = 0; j < simplexSize; ++j)
					{
						pts[j] = V[F[i][j]];
					}
					curA = tet_signed_volume(pts);
					minA = (curA < minA) ? curA : minA;
				}
				//tet_signed_volume returns 6 times of the volume
				minA = minA/6;
			}
			minAreaRecord.push_back(minA);
		}
	}

	//custom stop criteria
	bool stopQ()
	{
		if (stopCode == "all_good") return stopQ_all_good();
		//default
		return false;
	}

	bool stopQ_all_good()
	{
		bool good = true;
		unsigned nF = F.size();
		unsigned simplexSize = F[0].size();
		std::vector<std::vector<double> > pts(simplexSize);
		if (simplexSize == 3)
		{ //triMesh
			for (int i = 0; i < nF; ++i)
			{
				for (int j = 0; j < simplexSize; ++j)
				{
					pts[j] = V[F[i][j]];
				}
				if (tri_signed_area(pts) <= 0)
				{
					good = false;
					break;
				}
			}
		}
		else { //tetMesh
			for (int i = 0; i < nF; ++i)
			{
				for (int j = 0; j < simplexSize; ++j)
				{
					pts[j] = V[F[i][j]];
				}
				if (tet_signed_volume(pts) <= 0)
				{
					good = false;
					break;
				}
			}
		}
		
		return good;
	}

	//export result
	bool exportResult(const char* filename)
	{
		std::ofstream out_file(filename);
		if (! out_file.is_open()) {
			std::cerr << "Failed to open " << filename << "!" << std::endl;
			return false;
		}

		//precision of output
		typedef std::numeric_limits< double > dbl;
		out_file.precision(dbl::max_digits10);

    	//write the file
		unsigned nv = V.size();
		unsigned ndim = V[0].size();
		out_file << "resV " << nv << " " << ndim << "\n";
		for (int i = 0; i < nv; ++i)
		{
			for (int j = 0; j < ndim; ++j)
			{
				out_file << V[i][j] << " ";
			}
		}
		out_file << std::endl;


		if (record_vert)
		{
			unsigned n_record = vertRecord.size();
			nv = vertRecord[0].size();
			ndim = vertRecord[0][0].size();
			out_file << "vert " << n_record << " " << nv << " " << ndim << std::endl;
			for (int i = 0; i < n_record; ++i)
			{
				for (int j = 0; j < nv; ++j)
				{
					for (int k = 0; k < ndim; ++k)
					{
						out_file << vertRecord[i][j][k] << " ";
					}
				}
			}
			out_file << std::endl;
		}

		if (record_energy)
		{
			unsigned n_record = energyRecord.size();
			out_file << "energy " << n_record << std::endl;
			for (int i = 0; i < n_record; ++i)
			{
				out_file << energyRecord[i] << " ";
			}
			out_file << std::endl;
		}
		
		if (record_minArea)
		{
			unsigned n_record = minAreaRecord.size();
			out_file << "minArea " << n_record << std::endl;
			for (int i = 0; i < n_record; ++i)
			{
				out_file << minAreaRecord[i] << " ";
			}
			out_file << std::endl;
		}

		out_file.close();
		return true;
	}
};


//
double lifted_func(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
	LiftedData *data = (LiftedData *) my_func_data;

	if (data->solutionFound)
	{
		if (!grad.empty())
		{
			for (int i = 0; i < grad.size(); ++i)
			{
				grad[i] = 0;
			}
		}
		return data->lastFunctionValue;
	}

	std::vector<std::vector<double> > &V = data->V;  // nV x ndim
	unsigned nV = V.size();
	unsigned ndim = V[0].size();
	const std::vector<unsigned> &freeI = data->freeI;  //free vert indices

	const std::vector<std::vector<unsigned> > &F = data->F; // nF x simplexSize
	unsigned nF = F.size();
	unsigned simplexSize = F[0].size();
	unsigned nEdge = (simplexSize*(simplexSize-1))/2;

	//update V with x
	for (int i = 0; i < freeI.size(); ++i)
	{
		for (int j = 0; j < ndim; ++j)
		{
			V[freeI[i]][j] = x[ndim*i+j];
		}
	}

	// custom stop criterion
	if (data->stopQ())
	{
		data->solutionFound = true;
	}

	//compute squared edge length
	std::vector<std::vector<double> > D;

	for (int i = 0; i < nF; ++i)
	{
		std::vector<double> Di;
		Di.reserve(nEdge);
		for (int j = 0; j < simplexSize; ++j)
		{
			for (int k = j+1; k < simplexSize; ++k)
			{
				Di.push_back(distance2(V[F[i][j]],V[F[i][k]]));
			}
		}
		D.push_back(Di);
	}

	mat_add_to(D, data->restD);

	//compute total unsigned volume
	std::vector<double> A(nF);
	double f = 0;
	if (simplexSize == 3)
	{ //triMesh
		for (int i = 0; i < nF; ++i)
		{
			A[i] = tri_area(D[i]);
			f += A[i];
		}
		f = f/4;
	}
	else { //tetmesh
		for (int i = 0; i < nF; ++i)
		{
			A[i] = tet_volume(D[i]);
			f += A[i];
		}
		f = f/12;
	}
	
	//test
	data->nb_feval += 1;

	//compute gradient
	if (!grad.empty())
	{
		//
		std::vector<std::vector<double> > dAdD;
		if (simplexSize == 3)
		{ //triMesh
			for (int i = 0; i < nF; ++i)
			{
				std::vector<double> gi = tri_grad(D[i]);
				double s = 1/(4*A[i]);
				vec_scale(gi,s);
				dAdD.push_back(gi);
			}
		}
		else { //tetMesh
			for (int i = 0; i < nF; ++i)
			{
				std::vector<double> gi = tet_grad(D[i]);
				double s = 1/(12*A[i]);
				vec_scale(gi,s);
				dAdD.push_back(gi);
			}
		}

		// init Grad matrix
		std::vector<std::vector<double> > G(nV);
		for (int i = 0; i < nV; ++i)
		{
			G[i].resize(ndim,0);
		}

		// fill Grad
		for (int i = 0; i < nF; ++i)
		{
			unsigned edge_idx = 0;
			for (int j = 0; j < simplexSize; ++j)
			{
				for (int k = j+1; k < simplexSize; ++k)
				{
					unsigned idx_j = F[i][j];
					unsigned idx_k = F[i][k];
					const std::vector<double> &vj = V[idx_j];
					const std::vector<double> &vk = V[idx_k];
					std::vector<double> ejk = vec_substract(vj,vk);
					double s = dAdD[i][edge_idx];
					vec_scale_add_to(G[idx_j],ejk,s);
					vec_scale_add_to(G[idx_k],ejk,-s);
					++edge_idx;
				}
			}
		}
		
		// grad = G[freeI]
		for (int i = 0; i < freeI.size(); ++i)
		{
			const std::vector<double> &gi = G[freeI[i]];
			for (int j = 0; j < ndim; ++j)
			{
				grad[ndim*i+j] = gi[j];
			}
		}

		//test
		data->nb_geval += 1;
	}

	//record information
	data->lastFunctionValue = f;
	data->record();
	return f;
}


//
int main(int argc, char const *argv[])
{
	const char* dataFile = (argc > 1) ? argv[1] : "./input";
	const char* optFile  = (argc > 2) ? argv[2] : "./solver_options";
	const char* resFile  = (argc > 3) ? argv[3] : "./result";

	//import data
	std::vector<std::vector<double> > restV;
	std::vector<std::vector<double> > initV;
	std::vector<std::vector<unsigned> > F;
	std::vector<unsigned> handles;

    bool succeed = importData(dataFile,restV,initV,F,handles);
    if (!succeed) {
        std::cout << "usage: findInjective [inputFile] [solverOptionsFile] [resultFile]" << std::endl;
        return -1;
    }

    //import options
    NloptOptionManager options(optFile);
    //std::cout << "--- options ---" << std::endl;
    //options.printOptions();


	//init
    LiftedData data(restV,initV,F,handles,options.form,options.alphaRatio,options.alpha);

	unsigned nv = restV.size();
	unsigned nfree = nv - handles.size();
	unsigned embed_dim = initV[0].size();
	unsigned problem_dim = nfree * embed_dim;

	//set algorithm
	nlopt::algorithm algo = nlopt::LD_LBFGS;
	if (options.algorithm == "LD_LBFGS")
	{
		algo = nlopt::LD_LBFGS;
	}

	nlopt::opt opt(algo, problem_dim);
	
	//set stop criteria
	opt.set_ftol_abs(options.ftol_abs);
	opt.set_ftol_rel(options.ftol_rel);
	opt.set_xtol_abs(options.xtol_abs);
	opt.set_xtol_rel(options.xtol_rel);
	opt.set_maxeval(options.maxeval);

	//pass relevant options to LiftedData
	data.stopCode = options.stopCode;
	data.set_record_flags(options.record);

	//
	opt.set_min_objective(lifted_func, &data);

	std::vector<double> x = data.x0;
	double minf;

	//optimize
	try{
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		nlopt::result result = opt.optimize(x, minf);
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Time difference: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " [microseconds]" << std::endl;
		std::cout << "result: ";
		switch(result) {
			case nlopt::SUCCESS:
				std::cout << "SUCCESS" << std::endl;
				break;
			case nlopt::STOPVAL_REACHED:
				std::cout << "STOPVAL_REACHED" << std::endl;
				break;
			case nlopt::FTOL_REACHED:
				std::cout << "FTOL_REACHED" << std::endl;
				break;
			case nlopt::XTOL_REACHED:
				std::cout << "XTOL_REACHED" << std::endl;
				break;
			case nlopt::MAXEVAL_REACHED:
				std::cout << "MAXEVAL_REACHED" << std::endl;
				break;
			case nlopt::MAXTIME_REACHED:
				std::cout << "MAXTIME_REACHED" << std::endl;
				break;
			default:
				std::cout << "unexpected return code!" << std::endl;
				break;
		}
		std::cout << "met custom stop criteria (" << options.stopCode << "): ";
		if (data.solutionFound) std::cout << "yes" << std::endl;
		else std::cout << "no" << std::endl;
		//
		//std::cout << data.nb_feval << " function evalations, ";
		//std::cout << data.nb_geval << " gradient evalations." << std::endl;

	}
	catch(std::exception &e) {
		std::cout << "nlopt failed: " << e.what() << std::endl;
	}

	//export result
	data.exportResult(resFile);

	return 0;
}
