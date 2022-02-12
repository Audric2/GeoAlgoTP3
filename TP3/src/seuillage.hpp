#ifndef SEUILLAGE_HPP
#define SEUILLAGE_HPP

#include "for_CGAL.hpp"
#include "otsu.hpp"

// Effectue un seuillage selon la moyenne de la propriete
Facet_int_map seuillageSimple(Polyhedron & mesh, Facet_double_map & propriete){
	Facet_int_map segmentation;
	double sumMoy = 0;
	int nbFaces = 0;
	// Calcul de la moyenne
	for (  Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) {
		sumMoy+= propriete[i];
		nbFaces++;
	}
	sumMoy/=nbFaces;
	// Seuillage en utilisant la moyenne comme seuil
	for (  Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		segmentation[face] = (propriete[face]>sumMoy)?0:1;
	}
	return segmentation;
}

// Effectue un seuillage en coupant en n selon l'ordre des propriete
Facet_int_map seuillageMultiple(Polyhedron & mesh, Facet_double_map & propriete,int n){
	Facet_int_map segmentation;
	int nbFaces = 0;
	// Creation d'un vecteur avec toutes les valeurs triees dans l'ordre croissant
	std::vector<double> vValeurs(0);
	for (  Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) {
		vValeurs.push_back(propriete[i]);
		nbFaces++;	
	}
	std::sort (vValeurs.begin(), vValeurs.end());
	// Creation du vecteur de seuils  
	std::vector<double> seuils(n-1);
	for (int i = 0;i<n-1;i++) {
		int index = (i+1) * (nbFaces-1) / n;
		seuils[i] = vValeurs[index];
	}
	// Ne gere pas les egalite de seuils

	for (  Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		int index = 0;
		// On cherche la valeur la plus proche
		while(index < n-1 && propriete[face]>seuils[index]){// Dichotomie
			index++;
		}
		// On donne une classe entre 0 et n-1 en fonction de la position de la propriete dans le tableau
		segmentation[face] = index;
	}
	return segmentation;
}

// Effectue une segmentation en suivant la methode d'otsu
Facet_int_map seuillageOtsu(Polyhedron & mesh, Facet_double_map & propriete,int n = 64){
	Facet_int_map segmentation;
	// Creation d'un histogramme de taille n
	std::vector<int> histogramme(n);
	Facet_iterator face = mesh.facets_begin();
	double maxHisto = propriete[face], minHisto = propriete[face],rangeHisto = 1;
	// On cherche le minimum et le maximum de l'histo
	for (++face  ; face != mesh.facets_end(); ++face) {
		maxHisto = std::max(propriete[face], maxHisto);
		minHisto = std::min(propriete[face], minHisto);
	}
	//Verification que toutes les propriete ne valent pas la meme chose
	if(minHisto != maxHisto){
		rangeHisto = maxHisto - minHisto;
	}
	// On remplie l'histogramme grace au propriete des faces
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		int index = round((double)(n-1) * (propriete[face] - minHisto ) / rangeHisto);
		histogramme[index]++;
	}
	// On recupere lindex du seuil grace a la methode d'otsu
	int thresholdindex = otsu(histogramme);
	// On convertit l'index de seuil en une valeur de seuil 
	double threshold = rangeHisto * (double)(thresholdindex+0.5) / (double)(n-1) + minHisto;
	// Seuillage en utilisant la valeur donnee par otsu comme seuil
	for (  Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		segmentation[face] = propriete[face]>threshold;
	}
	return segmentation;
}

#endif

