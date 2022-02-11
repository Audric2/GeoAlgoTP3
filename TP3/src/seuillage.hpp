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
Facet_int_map seuillageMultiple(Polyhedron & mesh, Facet_double_map & propiete,int n){
	Facet_int_map segmentation;
	int nbFaces = 0;
	// Creation d'un vecteur avec toutes les valeurs triees dans l'ordre croissant
	std::vector<double> vValeurs(0);
	for (  Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) {
		vValeurs.push_back(propiete[i]);
		nbFaces++;	
	}
	std::sort (vValeurs.begin(), vValeurs.end());  
	
	for (  Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		int indexB = 0,indexHaut = nbFaces-1,indexMid;
		// On cherche la valeur la plus proche
		while(indexB+1<indexHaut){// Dichotomie
			indexMid = (indexB+indexHaut)/2;
			if(propiete[face]==vValeurs[indexMid]){
				indexB = indexMid;
				indexHaut = indexMid;
			}else if(propiete[face]<vValeurs[indexMid]){
				indexHaut = indexMid-1;
			}else{
				indexB = indexMid+1;
			}
		}
		// On donne une classe entre 0 et n-1 en fonction de la position de la propiete dans le tableau
		segmentation[face] = round(indexB * (n-1) / (nbFaces-1));
	}
	return segmentation;
}

// Effectue une segmentation en suivant la methode d'otsu
Facet_int_map seuillageOtsu(Polyhedron & mesh, Facet_double_map & propiete,int n = 64){
	Facet_int_map segmentation;
	// Creation d'un histogramme de taille n
	std::vector<int> histogramme(n);
	Facet_iterator face = mesh.facets_begin();
	double maxHisto = propiete[face], minHisto = propiete[face],rangeHisto = 1;
	// On cherche le minimum et le maximum de l'histo
	for (++face  ; face != mesh.facets_end(); ++face) {
		maxHisto = std::max(propiete[face], maxHisto);
		minHisto = std::min(propiete[face], minHisto);
	}
	//Verification que toutes les propriete ne valent pas la meme chose
	if(minHisto != maxHisto){
		rangeHisto = maxHisto - minHisto;
	}
	// On remplie l'histogramme grace au propriete des faces
	for (Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		int index = round((double)(n-1) * (propiete[face] - minHisto ) / rangeHisto);
		histogramme[index]++;
	}
	// On recupere lindex du seuil grace a la methode d'otsu
	int thresholdindex = otsu(histogramme);
	// On convertit l'index de seuil en une valeur de seuil 
	double threshold = rangeHisto * (double)(thresholdindex+0.5) / (double)(n-1) + minHisto;
	// Seuillage en utilisant la valeur donnee par otsu comme seuil
	for (  Facet_iterator face = mesh.facets_begin(); face != mesh.facets_end(); ++face) {
		segmentation[face] = propiete[face]>threshold;
	}
	return segmentation;
}

#endif

