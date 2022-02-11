#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/squared_distance_3.h> 
#include <iostream>
#include <fstream>
#include <random>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
typedef Polyhedron::Facet_handle Facet_handle;
typedef std::map<Facet_handle, double> Facet_double_map;
typedef std::map<Facet_handle, int> Facet_int_map;

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
	int n = 4;
	std::vector<int> histogramme(n);
	Facet_iterator i = P.facets_begin();
	double maxHisto = valeurs[i];
	double minHisto = valeurs[i];
	for (++i  ; i != P.facets_end(); ++i) {
		maxHisto = std::max(valeurs[i], maxHisto);
		minHisto = std::min(valeurs[i], minHisto);
	}
	for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		int index = round((double)(n-1) * (valeurs[i] - minHisto ) / (maxHisto - minHisto));
		histogramme[index]++;
	}
	for(int index = 0;index<n;index++){
		printf("%d,",histogramme[index]);
	}
	printf("\n");
	/* a faire 
	for(int index = 0;index<n;index++){
		int w1=0,w2=0;
		for(int index = 0;index<n;index++){

		}
	}*/
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

	std::ofstream output(std::string(argv[1]) + "colored.off");
	Facet_double_map per = aire(mesh);
	//Facet_int_map segmentation = seuillageSimple(mesh, per);
	Facet_int_map segmentation = seuillageOtsu(mesh, per);
	Facet_int_map classes = segmentationParCC(mesh, segmentation);
	saveClasses(output, mesh, classes);
	output.close();
	return 0;
}
