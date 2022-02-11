#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/squared_distance_3.h> 
#include <iostream>
#include <fstream>
#include <random>
#include "otsu_histogramme.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
typedef Polyhedron::Facet_handle Facet_handle;
typedef std::map<Facet_handle, double> Facet_double_map;
typedef std::map<Facet_handle, int> Facet_int_map;


double min(double a, double b){
    return (a<b)?a:b;
}

double max(double a, double b){
    return (a>b)?a:b;
}

void savePer(std::ofstream & file, Polyhedron & P, Facet_double_map & perimetres){
	// Write polyhedron in Object File Format (OFF).
    CGAL::set_ascii_mode(file);
    file << "OFF" << std::endl << P.size_of_vertices() << ' '
              << P.size_of_facets() << " 0" << std::endl;
    std::copy( P.points_begin(), P.points_end(),
               std::ostream_iterator<Point_3>(file, "\n"));

	double maxPerimetre = 0;
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		maxPerimetre = std::max(perimetres[i],maxPerimetre);
	}
    for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_around_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        file << CGAL::circulator_size(j) << ' ';
        do {
            file << ' ' << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
		file << " " << perimetres[i]/maxPerimetre << " " << 0 << " " << 0 << " " << 0.75;
        file << std::endl;
    }
}

void saveClasses(std::ofstream & file, Polyhedron & P, Facet_int_map & classes){
	// Write polyhedron in Object File Format (OFF).
    CGAL::set_ascii_mode(file);
    file << "OFF" << std::endl << P.size_of_vertices() << ' '
              << P.size_of_facets() << " 0" << std::endl;
    std::copy( P.points_begin(), P.points_end(),
               std::ostream_iterator<Point_3>(file, "\n"));

	int nbClasses = 0;
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		nbClasses = std::max(classes[i],nbClasses);
	}
	nbClasses++;
	std::vector<int> vecCol = std::vector<int>(nbClasses);

	std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist24(0,0xFFFFFF);
	dist24(rng);
	for(int i = 0;i<nbClasses;++i){
		vecCol[i] = dist24(rng); 
	}

    for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_around_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        file << CGAL::circulator_size(j) << ' ';
        do {
            file << ' ' << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
		//std::cout << vecCol[classes[i]] << "\n";
		file << " " << ((vecCol[classes[i]]>>16)&0xFF)/255.f << " " << ((vecCol[classes[i]]>>8)&0xFF)/255.f << " " << ((vecCol[classes[i]])&0xFF)/255.f << " " << 0.75;
        file << std::endl;
    }
}

Facet_int_map seuillageSimple(Polyhedron & P, Facet_double_map & perimetre){
	Facet_int_map classes;
	double Moy = 0;
	int nbFaces = 0;
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		Moy+= perimetre[i];
		nbFaces++;
    }
	Moy/=nbFaces;
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		classes[i] = (perimetre[i]>Moy)?0:1;
    }
	return classes;
}

Facet_int_map seuillageMultiple(Polyhedron & P, Facet_double_map & valeurs,int n){
	Facet_int_map classes;
	int nbFaces = 0;
	std::vector<double> vValeurs(0);
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		vValeurs.push_back(valeurs[i]);
	    nbFaces++;	
    }
    
	std::sort (vValeurs.begin(), vValeurs.end());  
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		int indexB = 0,indexHaut = nbFaces-1,indexMid;
		while(indexB+1!=indexHaut){
			indexMid = (indexB+indexHaut)/2;
			if(valeurs[i]==vValeurs[indexMid]){
				indexB = indexMid;
				indexHaut = indexMid+1;
			}else if(valeurs[i]<vValeurs[indexMid]){
				indexHaut = indexMid-1;
			}else{
				indexB = indexMid+1;
			}
		}
		classes[i] = round(indexB * (n-1) / (nbFaces-1));
    }
	return classes;
}

void parcoursFaces(int classe, Facet_handle & i, Facet_int_map & classes,Facet_int_map & segmentation){
	auto result = classes.find(i);
    if (result != classes.end()) {
		return;
    }
	classes[i] = classe;
	Halfedge_around_facet_circulator heIt = i->facet_begin();
	do{
		auto faceAdj = heIt->opposite()->facet();
		if(faceAdj!=NULL && segmentation[i] == segmentation[faceAdj]){
			parcoursFaces(classe, faceAdj, classes, segmentation);
		}
		++heIt;
	}while(heIt != i->facet_begin());
}

Facet_int_map segmentationParCC(Polyhedron & P, Facet_int_map & segmentation){
	Facet_int_map classes;
	int classe = 0;

	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		parcoursFaces(classe++, i, classes, segmentation);
    }
	return classes;
}

Facet_int_map seuillageOtsu(Polyhedron & P, Facet_double_map & valeurs){
	Facet_int_map segmentation;
	int n = 64;
	std::vector<int> histogramme(n);
	Facet_iterator i = P.facets_begin();
	double maxHisto = valeurs[i], minHisto = valeurs[i],L = 1;
	for (++i  ; i != P.facets_end(); ++i) {
		maxHisto = std::max(valeurs[i], maxHisto);
		minHisto = std::min(valeurs[i], minHisto);
	}
	if(minHisto != maxHisto){
	    L = maxHisto - minHisto;
	}
	for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		int index = round((double)(n-1) * (valeurs[i] - minHisto ) / L);
		histogramme[index]++;
	}
	/*
	for(int index = 0;index<n;index++){
		printf("%d,",histogramme[index]);
	}*/
	
	int thresholdindex = otsu(histogramme);
	double threshold = L * (double)(thresholdindex+0.5) / (double)(n-1) + minHisto;
	int nb0 = 0,nb1 = 0;
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		segmentation[i] = valeurs[i]>threshold;
		if(segmentation[i]) nb1++;
		else nb0++;
    }
    printf("\n%d classe 0, %d classe 1\n",nb0,nb1);
	return segmentation;
}

std::map<Facet_handle, double> perimetre(Polyhedron & mesh){
	std::map<Facet_handle, double> perimetres;
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		double perimetre = 0;

		Halfedge_around_facet_circulator j = face->facet_begin();

		do{
			perimetre += CGAL::squared_distance(j->vertex()->point(),j->opposite()->vertex()->point());
		} while(++j != face->facet_begin());
		perimetres[face] = perimetre;
	}
	return perimetres;
}

double aireFace(Facet_handle face){
	double area = 0;
	Halfedge_around_facet_circulator j = face->facet_begin();
	auto v0 = j->vertex()->point();
	j++;
	auto v1 = j->vertex()->point();
	auto v2 = j->vertex()->point();
	j++;
	do{
		v1 = v2;
		v2 = j->vertex()->point();
		j++;
		area += sqrt(CGAL::squared_area(v0,v1,v2));
	} while(j != face->facet_begin());
	return area;
}

std::map<Facet_handle, double> aire(Polyhedron & mesh){
	std::map<Facet_handle, double> aires;
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		aires[face] = aireFace(face);
	}
	return aires;
}

double angleMinFace(Facet_handle face){
	double minAngle = 360;
	Halfedge_around_facet_circulator j = face->facet_begin();
	do{
		auto lj = j;
		j++;
		double angle =CGAL:: approximate_angle(
		    Vector_3(lj->opposite()->vertex()->point(),
		    lj->vertex()->point()),
		    Vector_3 (j->vertex()->point(),
		    j->opposite()->vertex()->point()));
	    if(angle!=0){
	        //printf("%lf,",angle);
	        if(angle<0){
	            angle = abs(angle);
	        }
	        if(angle>180){
	            angle = 360-angle;
	        }
	        minAngle = min(minAngle,angle);
	    }
	} while(j != face->facet_begin());
	return minAngle;
}

std::map<Facet_handle, double> minAngle(Polyhedron & mesh){

	std::map<Facet_handle, double> angles;
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		angles[face] = angleMinFace(face);
		if(angles[face]>60){
		    printf("ATTENTION LE MESH N'EST PAS COMPOSE QUE DE TRIANGLE\n");
		    printf("%.2lf\n",angles[face]);
		}
	}
	return angles;
}

std::map<Facet_handle, double> angleFaceV(Polyhedron & mesh,Vector_3 & v){

	std::map<Facet_handle, double> angles;
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
	    Halfedge_around_facet_circulator j = face->facet_begin();
		auto lj = j;
		j++;
		auto normal = CGAL::normal(
		    lj->opposite()->vertex()->point(),
		    lj->vertex()->point(),
		    j->vertex()->point());
        double angle = CGAL::approximate_angle(normal, v);
	    angles[face] = angle;
	}
	return angles;
}

std::map<Facet_handle, double> moyAngleVoisin(Polyhedron & mesh){
	std::map<Facet_handle, double> angles;
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
	    double angleMoy = 0;
	    int nbVoisin = 0;
	    Halfedge_around_facet_circulator j = face->facet_begin();
		auto lj = j;
		j++;
		auto normal = CGAL::normal(
		    lj->opposite()->vertex()->point(),
		    lj->vertex()->point(),
		    j->vertex()->point());
        Halfedge_around_facet_circulator heIt = face->facet_begin();
	    do{
		    auto faceAdj = heIt->opposite()->facet();
		    if(faceAdj!=NULL){
			    Halfedge_around_facet_circulator jF = faceAdj->facet_begin();
		        auto ljF = jF;
		        jF++;
		        auto normalF = CGAL::normal(
		            ljF->opposite()->vertex()->point(),
		            ljF->vertex()->point(),
		            jF->vertex()->point());
		        double angle = CGAL::approximate_angle(normal, normalF);
		        angleMoy += angle;
		        nbVoisin++;
		    }
		    ++heIt;
	    }while(heIt != face->facet_begin());
	    angles[face] = angleMoy/nbVoisin;
	}
	return angles;
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
		return 1;
	}
	
	Polyhedron mesh;
	std::ifstream input(argv[1]);
		if (!input || !(input >> mesh) || mesh.is_empty()) {
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}
  
  	unsigned int nbVerts = 0;
	for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
		++nbVerts;
	}
	std::cout << "Nombre de sommets: " << nbVerts << std::endl;
	
	unsigned int nbEdges = 0;
	for (Halfedge_iterator i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i) {
		++nbEdges;
	}
	nbEdges /= 2;
	std::cout << "Nombre d'arêtes: " << nbEdges << std::endl;

	unsigned int nbFaces = 0;
	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) {
		++nbFaces;
	}
	std::cout << "Nombre de faces: " << nbFaces << std::endl;
	
	unsigned int euler = nbVerts - nbEdges + nbFaces;
	unsigned int genus = (2 - euler) / 2;
	std::cout << "En supposant que le maillage contienne une unique surface sans bord, alors son genre est de " << genus << std::endl;
    Vector_3 v(1,0,0);
	Facet_double_map per = angleFaceV(mesh,v);
	Facet_int_map segmentation = seuillageMultiple(mesh, per, 4);//(mesh, per);
	Facet_int_map classes = segmentationParCC(mesh, segmentation);
	
	std::ofstream output(std::string(argv[1]) + "colored.off");
	saveClasses(output, mesh, classes);
	output.close();
	return 0;
}
