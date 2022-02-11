#ifndef SAVE_TO_FILE_HPP
#define SAVE_TO_FILE_HPP

/**
 * Fichier ou sont defini les sauveagardes dans des fichiers
 */

#include "for_CGAL.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** 
 * Sauvegarde le Polyhedron au format OFF 
 * donne une couleur a chaque face en fonction de son perimetre
*/
void savePer(std::ofstream & file, Polyhedron & P, Facet_double_map & perimetres){
	// Write polyhedron in Object File Format (OFF).
	CGAL::set_ascii_mode(file);
	// Entete du fichier .OFF
	file << "OFF" << std::endl << P.size_of_vertices() << ' '
			  << P.size_of_facets() << " 0" << std::endl;
	// Copie de toutes les vertices dans le fichier
	std::copy( P.points_begin(), P.points_end(),
			   std::ostream_iterator<Point_3>(file, "\n"));

	// Trouve le perimetre maximum
	double maxPerimetre = 0;
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		maxPerimetre = std::max(perimetres[i],maxPerimetre);
	}

	// Boucle sur toutes les faces
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		Halfedge_around_facet_circulator j = i->facet_begin();
		// Facets in polyhedral surfaces are at least triangles.
		CGAL_assertion( CGAL::circulator_size(j) >= 3);
		// Ecriture du nombre de points de la face
		file << CGAL::circulator_size(j) << ' ';
		do {
			 // Ecriture du numero de chaque point de la face
			file << ' ' << std::distance(P.vertices_begin(), j->vertex());
		} while ( ++j != i->facet_begin());
		// Ecriture d'une couleur RGB dependant du perimetre
		file << " " << perimetres[i]/maxPerimetre << " " << 0 << " " << 0 << " " << 0.75;
		file << std::endl;
	}
}

/** 
 * Sauvegarde le Polyhedron au format OFF 
 * donne une couleur a chaque face en fonction de sa classe
*/
void saveClasses(std::ofstream & file, Polyhedron & P, Facet_int_map & classes){
	// Write polyhedron in Object File Format (OFF).
	CGAL::set_ascii_mode(file);
	// Entete du fichier .OFF
	file << "OFF" << std::endl << P.size_of_vertices() << ' '
			  << P.size_of_facets() << " 0" << std::endl;
	// Copie de toutes les vertices dans le fichier
	std::copy( P.points_begin(), P.points_end(),
			   std::ostream_iterator<Point_3>(file, "\n"));
	
	// On regarde le nombre de classes dans le map 
	int nbClasses = 0;
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		nbClasses = std::max(classes[i],nbClasses);
	}
	nbClasses++;
	// Initialisation d'un vetceur de int pour stocker les couleurs
	std::vector<int> vecCol = std::vector<int>(nbClasses);

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> dist24(0,0xFFFFFF);
	// Initialisation d'une couleur aleatoire par classe
	for(int i = 0;i<nbClasses;++i){
		vecCol[i] = dist24(rng); 
	}
	
	// Boucle sur toutes les faces
	for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		Halfedge_around_facet_circulator j = i->facet_begin();
		// Facets in polyhedral surfaces are at least triangles.
		CGAL_assertion( CGAL::circulator_size(j) >= 3);
		// Ecriture du nombre de points de la face
		file << CGAL::circulator_size(j) << ' ';
		do {
			// Ecriture du numero de chaque point de la face
			file << ' ' << std::distance(P.vertices_begin(), j->vertex());
		} while ( ++j != i->facet_begin());
		// Ecriture d'une couleur RGB dependant de la classe
		file << " " << ((vecCol[classes[i]]>>16)&0xFF)/255.f << " " << ((vecCol[classes[i]]>>8)&0xFF)/255.f << " " << ((vecCol[classes[i]])&0xFF)/255.f << " " << 0.75;
		file << std::endl;
	}
}

#endif

