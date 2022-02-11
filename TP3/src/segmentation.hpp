#ifndef SEGMENTATION_HPP
#define SEGMENTATION_HPP

#include "for_CGAL.hpp"

// Rassemblement des faces connexes grace a un parcours des faces voisines recursivement 
void parcoursFaces(int classe, Facet_handle & face, Facet_int_map & classes,Facet_int_map & segmentation){
	// On verifie si la face a deja une classe
	auto result = classes.find(face);
    if (result != classes.end()) {
		return;
    }
	// On affecte une classe a la face
	classes[face] = classe;
	// Puis on relance la fonction depuis tous les voisins qui ont la meme segmentation que nous
	Halfedge_around_facet_circulator halfedgeIt = face->facet_begin();
	do{
		// Recuperation d'une face adjacente
		auto faceAdjacente = halfedgeIt->opposite()->facet();
		// Verification que la face existe et quelle a la meme segmentation que nous 
		if(faceAdjacente!=NULL && segmentation[face] == segmentation[faceAdjacente]){
			// Relance le parcours depuis cette face
			parcoursFaces(classe, faceAdjacente, classes, segmentation);
		}
	}while(++halfedgeIt != face->facet_begin());
}

// Regroupement des faces connexes quand elle sont de la meme segmentation
Facet_int_map segmentationParCC(Polyhedron & mesh, Facet_int_map & segmentation){
	Facet_int_map classes;
	int classe = 0;
	// On lance un parcours depuis toutes les faces avec une classe qui augmente
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		parcoursFaces(classe++, face, classes, segmentation);
    }
	return classes;
}

#endif

