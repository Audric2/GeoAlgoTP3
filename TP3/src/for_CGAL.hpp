#ifndef FOR_CGAL_HPP
#define FOR_CGAL_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/squared_distance_3.h> 

// typedef utilises
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

#endif

