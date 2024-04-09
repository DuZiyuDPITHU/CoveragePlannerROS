/*
 * polygon_coverage_planning implements algorithms for coverage planning in
 * general polygons with holes. Copyright (C) 2019, Rik Bähnemann, Autonomous
 * Systems Lab, ETH Zürich
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "cgal_comm.h"
#include "visibility_graph.h"
#include "visibility_polygon.h"


namespace polygon_coverage_planning {
namespace visibility_graph {

VisibilityGraph::VisibilityGraph(const PolygonWithHoles& polygon)
    : GraphBase(), polygon_(polygon) {
  // Build visibility graph.
  is_created_ = create();
}

bool VisibilityGraph::create() {
  clear();
  // Sort vertices.
  sortVertices(&polygon_);
  // Select shortest path vertices.
  std::vector<VertexConstCirculator> graph_vertices;
  findConcaveOuterBoundaryVertices(&graph_vertices);
  findConvexHoleVertices(&graph_vertices);

  for (const VertexConstCirculator& v : graph_vertices) {
    // Compute visibility polygon.
    Polygon_2 visibility;
    if (!computeVisibilityPolygon(polygon_, *v, &visibility)) {
      std::cout<<"Cannot compute visibility polygon."<<std::endl;
      return false;
    }
    if (!addNode(NodeProperty(*v, visibility))) {
      return false;
    }
  }

//    std::cout<<"Created visibility graph with "
//                   << graph_.size() << " nodes and " << edge_properties_.size()
//                   << " edges."<<std::endl;

  return true;
}

void VisibilityGraph::findConcaveOuterBoundaryVertices(
    std::vector<VertexConstCirculator>* concave_vertices) const {

  VertexConstCirculator vit = polygon_.outer_boundary().vertices_circulator();
  do {
    Triangle_2 triangle(*std::prev(vit), *vit, *std::next(vit));
    CGAL::Orientation orientation = triangle.orientation();

    if (orientation == CGAL::CLOCKWISE) concave_vertices->push_back(vit);
  } while (++vit != polygon_.outer_boundary().vertices_circulator());
}

void VisibilityGraph::findConvexHoleVertices(
    std::vector<VertexConstCirculator>* convex_vertices) const {


  for (PolygonWithHoles::Hole_const_iterator hit = polygon_.holes_begin();
       hit != polygon_.holes_end(); ++hit) {

    VertexConstCirculator vit = hit->vertices_circulator();
    do {
      Triangle_2 triangle(*std::prev(vit), *vit, *std::next(vit));
      CGAL::Orientation orientation = triangle.orientation();

      if (orientation == CGAL::CLOCKWISE) convex_vertices->push_back(vit);
    } while (++vit != hit->vertices_circulator());
  }
}

bool VisibilityGraph::addEdges() {
  if (graph_.empty()) {
      std::cout<<"Cannot add edges to an empty graph."<<std::endl;
    return false;
  }

  const size_t new_id = graph_.size() - 1;
  for (size_t adj_id = 0; adj_id < new_id; ++adj_id) {
    const NodeProperty* new_node_property = getNodeProperty(new_id);
    const NodeProperty* adj_node_property = getNodeProperty(adj_id);
    if (adj_node_property == nullptr) {
        std::cout<<"Cannot access potential neighbor."<<std::endl;
      return false;
    }
    if (pointInPolygon(new_node_property->visibility,
                       adj_node_property->coordinates)) {
      EdgeId forwards_edge_id(new_id, adj_id);
      EdgeId backwards_edge_id(adj_id, new_id);
      const double cost = computeEuclideanSegmentCost(
          new_node_property->coordinates,
          adj_node_property->coordinates);  // Symmetric cost.
      if (!addEdge(forwards_edge_id, EdgeProperty(), cost) ||
          !addEdge(backwards_edge_id, EdgeProperty(), cost)) {
        return false;
      }
    }
  }
  return true;
}

bool VisibilityGraph::solve(const Point_2& start, const Point_2& goal,
                            std::vector<Point_2>* waypoints) const {

  waypoints->clear();

  // Make sure start and end are inside the polygon.
  const Point_2 start_new = pointInPolygon(polygon_, start)
                                ? start
                                : projectPointOnHull(polygon_, start);
  const Point_2 goal_new = pointInPolygon(polygon_, goal)
                               ? goal
                               : projectPointOnHull(polygon_, goal);

  // Compute start and goal visibility polygon.
  Polygon_2 start_visibility, goal_visibility;
  if (!computeVisibilityPolygon(polygon_, start_new, &start_visibility) ||
      !computeVisibilityPolygon(polygon_, goal_new, &goal_visibility)) {
    return false;
  }

  // Find shortest path.
  return solve(start_new, start_visibility, goal_new, goal_visibility,
               waypoints);
}

bool VisibilityGraph::solve(const Point_2& start,
                            const Polygon_2& start_visibility_polygon,
                            const Point_2& goal,
                            const Polygon_2& goal_visibility_polygon,
                            std::vector<Point_2>* waypoints) const {

  waypoints->clear();

  if (!is_created_) {
      std::cout<<"Visibility graph not initialized."<<std::endl;
    return false;
  } else if (!pointInPolygon(polygon_, start) ||
             !pointInPolygon(polygon_, goal)) {
      std::cout<<"Start or goal is not in polygon."<<std::endl;
    return false;
  }

  VisibilityGraph temp_visibility_graph = *this;
  // Add start and goal node.
  if (!temp_visibility_graph.addStartNode(
          NodeProperty(start, start_visibility_polygon)) ||
      !temp_visibility_graph.addGoalNode(
          NodeProperty(goal, goal_visibility_polygon))) {
    return false;
  }

  // Check if start and goal are in line of sight.
  size_t start_idx = temp_visibility_graph.getStartIdx();
  size_t goal_idx = temp_visibility_graph.getGoalIdx();
  const NodeProperty* start_node_property =
      temp_visibility_graph.getNodeProperty(start_idx);
  if (start_node_property == nullptr) {
    return false;
  }
  if (pointInPolygon(start_node_property->visibility, goal)) {
    waypoints->push_back(start);
    waypoints->push_back(goal);
    return true;
  }

  // Find shortest way using A*.
  Solution solution;
  if (!temp_visibility_graph.solveAStar(start_idx, goal_idx, &solution)) {
      std::cout<<
        "Could not find shortest path. Graph not fully connected."<<std::endl;
    return false;
  }

  // Reconstruct waypoints.
  return temp_visibility_graph.getWaypoints(solution, waypoints);
}

bool VisibilityGraph::getWaypoints(const Solution& solution,
                                   std::vector<Point_2>* waypoints) const {

  waypoints->resize(solution.size());
  for (size_t i = 0; i < solution.size(); i++) {
    const NodeProperty* node_property = getNodeProperty(solution[i]);
    if (node_property == nullptr) {
        std::cout<<"Cannot reconstruct solution."<<std::endl;
      return false;
    }
    (*waypoints)[i] = node_property->coordinates;
  }
  return true;
}

bool VisibilityGraph::calculateHeuristic(size_t goal,
                                         Heuristic* heuristic) const {

  heuristic->clear();

  const NodeProperty* goal_node_property = getNodeProperty(goal);
  if (goal_node_property == nullptr) {
      std::cout<<"Cannot find goal node property to calculate heuristic."<<std::endl;
    return false;
  }

  for (size_t adj_id = 0; adj_id < graph_.size(); ++adj_id) {
    const NodeProperty* adj_node_property = getNodeProperty(adj_id);
    if (adj_node_property == nullptr) {
        std::cout<<
          "Cannot access adjacent node property to calculate heuristic."<<std::endl;
      return false;
    }
    (*heuristic)[adj_id] = computeEuclideanSegmentCost(
        adj_node_property->coordinates, goal_node_property->coordinates);
  }

  return true;
}

bool VisibilityGraph::solveWithOutsideStartAndGoal(
    const Point_2& start, const Point_2& goal,
    std::vector<Point_2>* waypoints) const {


  if (solve(start, goal, waypoints)) {
    if (!pointInPolygon(polygon_, start)) {
      waypoints->insert(waypoints->begin(), start);
    }
    if (!pointInPolygon(polygon_, goal)) {
      waypoints->push_back(goal);
    }
    return true;
  } else {
    return false;
  }
}

double VisibilityGraph::computeEuclideanSegmentCost(const Point_2& from,
                                                    const Point_2& to) const {
  return std::sqrt(CGAL::to_double(Segment_2(from, to).squared_length()));
}

}  // namespace visibility_graph
}  // namespace polygon_coverage_planning
