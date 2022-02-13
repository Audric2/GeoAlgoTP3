#include "for_CGAL.hpp"

#include <iostream>
#include <fstream>
#include <random>
#include "save_to_file.hpp"
#include "seuillage.hpp"
#include "segmentation.hpp"
#include "local_property.hpp"

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
	for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) ++nbVerts;
	std::cout << "Nombre de sommets: " << nbVerts << std::endl;
	
	unsigned int nbEdges = 0;
	for (Halfedge_iterator i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i) ++nbEdges;
	nbEdges /= 2;
	std::cout << "Nombre d'arêtes: " << nbEdges << std::endl;
		
	unsigned int nbFaces = 0;
	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) ++nbFaces;
	std::cout << "Nombre de faces: " << nbFaces << std::endl;
	
	unsigned int euler = nbVerts - nbEdges + nbFaces;
	unsigned int genus = (2 - euler) / 2;
	std::cout << "En supposant que le maillage contienne une unique surface sans bord, alors son genre est de " << genus << std::endl;
	
	/*
		Test des fonctions implementees
	*/
	
	Vector_3 v(1,0,0);
	// Calcul d'une propriete locale pour chaque face
	Facet_double_map propriete = perimetre(mesh);
	// Seuillage de la dite propriete
	Facet_int_map segmentation = seuillageOtsu(mesh, propriete);
	// Creation d'une autre segmentation basee sur le seuillage precedent et les composantes connexes
	Facet_int_map classes = segmentationParCC(mesh, segmentation);

	// On sauvegarde le resultat dans un fichier .off	
	std::ofstream output(std::string(argv[1]) + "colored.off");
	saveClasses(output, mesh, classes);
	output.close();
	
	return 0;
}
