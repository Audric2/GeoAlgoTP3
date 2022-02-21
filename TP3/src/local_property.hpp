#ifndef LOCAL_PROPERTY_HPP
#define LOCAL_PROPERTY_HPP
/**
 * Fonction de calcul de proprietes locales sur les faces
*/
#include "for_CGAL.hpp"

// Calcule le perimetre de toutes les faces d'un Polyhedron
std::map<Facet_handle, double> perimetre(Polyhedron & mesh){
	// Creation de la map de stockage 
	std::map<Facet_handle, double> perimetres;
	// Boucle sur toutes les faces
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		double perimetre = 0;
		// Recuperation du circulateur de demi-arrete
		Halfedge_around_facet_circulator halfedgeIt = face->facet_begin();
		do{
			// Ajout de la longueur de l'arrete au perimetre
			perimetre += CGAL::squared_distance(halfedgeIt->vertex()->point(),halfedgeIt->opposite()->vertex()->point());
		} while(++halfedgeIt != face->facet_begin());
		// Enregistrement du perimetre associe a la face
		perimetres[face] = perimetre;
	}
	return perimetres;
}

// Calcul de l'aire d'une face
double aireFace(Facet_handle face){
	double area = 0;
	// Recuperation du circulateur de demi-arrete
	Halfedge_around_facet_circulator halfedgeIt = face->facet_begin();
	// Recuperation du premier vertex 
	Point_3 v0 = halfedgeIt->vertex()->point();halfedgeIt++;
	// Recuperation du deuxieme vertex 
	Point_3 v1,v2 = halfedgeIt->vertex()->point();halfedgeIt++;
	do{
		// On prend les points decrivant le prochain triangle
		v1 = v2;
		v2 = halfedgeIt->vertex()->point();
		// On calcule l'aire de ce triangle qu'on ajoute a l'aire total
		area += sqrt(CGAL::squared_area(v0,v1,v2));
	} while(++halfedgeIt != face->facet_begin());
	return area;
}

// Calcul l'aire de toutes les faces d'un Polyhedron
std::map<Facet_handle, double> aire(Polyhedron & mesh){
	std::map<Facet_handle, double> aires;
	// Boucle sur toutes les faces
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		// Association d'une face a son aire
		aires[face] = aireFace(face);
	}
	return aires;
}

// Calcul l'angle minimum d'une face
double angleMinFace(Facet_handle face){
	double minAngle = 360;
	// Recuperation du circulateur de demi-arrete
	Halfedge_around_facet_circulator halfedgeIt = face->facet_begin();
	do{
		// Recuperation de l'ancienne arrete et passage a la nouvelle
		auto lastHalfedgeIt = halfedgeIt++;
		// Calcul de l'angle entre l'ancienne arrete et la nouvelle
		double angle = CGAL:: approximate_angle(
			Vector_3(lastHalfedgeIt->opposite()->vertex()->point(),lastHalfedgeIt->vertex()->point()),
			Vector_3 (halfedgeIt->vertex()->point(),halfedgeIt->opposite()->vertex()->point()));
		if(angle){
			if(angle<0){// Si l'angle et negatif on le passe en positif
				angle = abs(angle);
			}
			// On verifie si l'angle est le nouveau minimum
			minAngle = std::min(minAngle,angle);
		}else printf("Cette face a un probleme, un des angle vaut 0\n");
	} while(halfedgeIt != face->facet_begin());
	return minAngle;
}

// Calcul l'angle minimum de toutes les faces d'un Polyhedron
std::map<Facet_handle, double> minAngle(Polyhedron & mesh){

	std::map<Facet_handle, double> angles;
	// Boucle sur toutes les faces
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		// Association d'une face a son angle minimum
		angles[face] = angleMinFace(face);
		if(angles[face]>60) printf("ATTENTION LE MESH N'EST PAS COMPOSE QUE DE TRIANGLE (Angle min : %.2lf > 60)\n",angles[face]);
	}
	return angles;
}

// Calcul l'angle entre la normale d'une face et un vecteur v pour toutes les faces d'un Polyhedron
std::map<Facet_handle, double> angleFaceV(Polyhedron & mesh,Vector_3 & v){

	std::map<Facet_handle, double> angles;
	// Boucle sur toutes les faces
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		// Recuperation du circulateur de demi-arrete
		Halfedge_around_facet_circulator halfedgeIt = face->facet_begin();
		auto lastHalfedgeIt = halfedgeIt++;
		// Calcul de la normal de la face
		auto normal = CGAL::normal(
			lastHalfedgeIt->opposite()->vertex()->point(),
			lastHalfedgeIt->vertex()->point(),
			halfedgeIt->vertex()->point());
		// Calcul de l'angle entre la normal et le vecteur v
		double angle = CGAL::approximate_angle(normal, v);
		// Association d'une face a son angle entre v et sa normal
		angles[face] = angle;
	}
	return angles;
}

// Calcul la moyenne des angles d'une face avec ses voisins pour toutes les faces d'un Polyhedron
std::map<Facet_handle, double> moyAngleVoisin(Polyhedron & mesh){
	std::map<Facet_handle, double> angles;
	// Boucle sur toutes les faces
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		double sumAngle = 0;
		int nbVoisin = 0;
		// Recuperation du circulateur de demi-arrete
		Halfedge_around_facet_circulator halfedgeIt = face->facet_begin();
		auto lastHalfedgeIt = halfedgeIt++;
		// Calcul de la normal de la face
		auto normal = CGAL::normal(
			lastHalfedgeIt->opposite()->vertex()->point(),
			lastHalfedgeIt->vertex()->point(),
			halfedgeIt->vertex()->point());

		// Remise a zero du circulateur de demi-arrete
		halfedgeIt = face->facet_begin();
		do{
			// Recuperation d'une face adjacente
			auto faceAdjacente = halfedgeIt->opposite()->facet();
			if(faceAdjacente!=NULL){ // Si on est pas sur un bord
				// Recuperation du circulateur de demi-arrete de la face adjacente
				Halfedge_around_facet_circulator halEdgeAdjacentIt = faceAdjacente->facet_begin();
				auto lasthalEdgeAdjacentIt = halEdgeAdjacentIt;
				halEdgeAdjacentIt++;
				// Calcul de la normale de la face adjacente
				auto normalF = CGAL::normal(
					lasthalEdgeAdjacentIt->opposite()->vertex()->point(),
					lasthalEdgeAdjacentIt->vertex()->point(),
					halEdgeAdjacentIt->vertex()->point());
				// Calcul de l'angle entre la normale et la normale de la face adjacente
				double angle = CGAL::approximate_angle(normal, normalF);
				// Ajout de l'angle pour la moyenne
				sumAngle += angle;
				nbVoisin++;
			}
		}while(++halfedgeIt != face->facet_begin());
		// Association d'une face a l'angle moyen entre sa normale et les normales des faces voisines
		angles[face] = sumAngle/nbVoisin;
	}
	return angles;
}


#endif

