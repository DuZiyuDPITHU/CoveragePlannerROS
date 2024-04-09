#include "coverage_planner.h"
#include <ros/ros.h>
#include <nav_msgs/Path.h>

std::vector<std::vector<Point_2>> polys;
Point_2 start(544, 447);
int contour_rcv_count = 0;

std::vector<Point_2> globalPathPlan(std::vector<std::vector<Point_2>> polys, Point_2 start)
{
    // compute main direction
    // Futher explaination: find on which direction the external contour length add up the most.

    // [0,180)
    std::vector<int> line_deg_histogram(180);
    double line_len; // weight
    double line_deg;
    int line_deg_idx;

    auto ext_poly = polys.front();
    ext_poly.emplace_back(ext_poly.front());
    for(int i = 1; i < ext_poly.size(); i++){
        int x_cur = CGAL::to_double(ext_poly[i].x());
        int x_pred= CGAL::to_double(ext_poly[i-1].x());
        int y_cur = CGAL::to_double(ext_poly[i].y());
        int y_pred= CGAL::to_double(ext_poly[i-1].y());
        line_len = std::sqrt(std::pow((x_cur-x_pred),2)+std::pow((y_cur-y_pred),2));
        // y-axis towards up, x-axis towards right, theta is from x-axis to y-axis
        line_deg = std::round(atan2(-(y_cur-y_pred),x_cur-x_pred)/M_PI*180.0); // atan2: (-180, 180]
        line_deg_idx = (int(line_deg) + 180) % 180; // [0, 180)
        line_deg_histogram[line_deg_idx] += int(line_len);
    }
    auto it = std::max_element(line_deg_histogram.begin(), line_deg_histogram.end());
    int main_deg = (it-line_deg_histogram.begin());
    std::cout<<"main deg: "<<main_deg<<std::endl;


    // construct polygon with holes
    ROS_INFO("global_planner: CONSTRUCTING POLYGON WITH HOLES.");
    std::vector<Point_2> outer_poly = polys.front();
    polys.erase(polys.begin());
    std::vector<std::vector<Point_2>> inner_polys = polys;

    Polygon_2 outer_polygon;
    printf("Outer polygon size: %ld\n", outer_poly.size());
    for(const auto& point : outer_poly){
        outer_polygon.push_back(point);
    }
    int num_holes = inner_polys.size();
    std::vector<Polygon_2> holes(num_holes);
    for(int i = 0; i < inner_polys.size(); i++){
        printf("Inner polygon %d size: %ld\n", i, inner_polys[i].size());
        for(const auto& point : inner_polys[i]){
            holes[i].push_back(point);
        }
    }
    PolygonWithHoles pwh(outer_polygon, holes.begin(), holes.end());
    // cell decomposition
    ROS_INFO("global_planner: CONDUCTING CELL DECOMPOSITION");
    std::vector<Polygon_2> bcd_cells;
    polygon_coverage_planning::computeBestBCDFromPolygonWithHoles(pwh, &bcd_cells);
    // construct adjacent graph
    auto cell_graph = calculateDecompositionAdjacency(bcd_cells);
    int starting_cell_idx = getCellIndexOfPoint(bcd_cells, start);
    // Futher explaination: conduct DFS to findout order of visiting cell nodes.
    ROS_INFO("global_planner: FINDING MOWER'S PATH");
    auto cell_idx_path = getTravellingPath(cell_graph, starting_cell_idx);
    std::cout<<"path length: "<<cell_idx_path.size()<<std::endl;
    std::cout<<"start";
    for(auto& cell_idx:cell_idx_path){
        std::cout<<"->"<<cell_idx;
    }
    std::cout<<std::endl;

    int sweep_step = 5;

    std::vector<std::vector<Point_2>> cells_sweeps;

    for (size_t i = 0; i < bcd_cells.size(); ++i) {
        // Compute all cluster sweeps.
        std::vector<Point_2> cell_sweep;
        Direction_2 best_dir;
        polygon_coverage_planning::findBestSweepDir(bcd_cells[i], &best_dir);
        polygon_coverage_planning::visibility_graph::VisibilityGraph vis_graph(bcd_cells[i]);

        bool counter_clockwise = true;
        polygon_coverage_planning::computeSweep(bcd_cells[i], vis_graph, sweep_step, best_dir, counter_clockwise, &cell_sweep);
        cells_sweeps.emplace_back(cell_sweep);
    }


    auto cell_intersections = calculateCellIntersections(bcd_cells, cell_graph);

    std::vector<Point_2> way_points;

    Point_2 point = start;
    std::list<Point_2> next_candidates;
    Point_2 next_point;
    std::vector<Point_2> shortest_path;

    if(doReverseNextSweep(start, cells_sweeps[cell_idx_path.front()])){
        shortest_path = getShortestPath(bcd_cells[cell_idx_path.front()], start, cells_sweeps[cell_idx_path.front()].back());
        way_points.insert(way_points.end(), shortest_path.begin(), std::prev(shortest_path.end()));
    } else{
        shortest_path = getShortestPath(bcd_cells[cell_idx_path.front()], start, cells_sweeps[cell_idx_path.front()].front());
        way_points.insert(way_points.end(), shortest_path.begin(), std::prev(shortest_path.end()));
    }

    point = way_points.back();

    for(size_t i = 0; i < cell_idx_path.size(); ++i){
        // has been cleaned?
        if(!cell_graph[cell_idx_path[i]].isCleaned){
            // need to reverse?
            if(doReverseNextSweep(point, cells_sweeps[cell_idx_path[i]])){
                way_points.insert(way_points.end(), cells_sweeps[cell_idx_path[i]].rbegin(), cells_sweeps[cell_idx_path[i]].rend());
            }else{
                way_points.insert(way_points.end(), cells_sweeps[cell_idx_path[i]].begin(), cells_sweeps[cell_idx_path[i]].end());
            }
            // now cleaned
            cell_graph[cell_idx_path[i]].isCleaned = true;
            // update current point
            point = way_points.back();
            // find shortest path to next cell
            if((i+1)<cell_idx_path.size()){
                next_candidates = cell_intersections[cell_idx_path[i]][cell_idx_path[i+1]];
                if(doReverseNextSweep(point, cells_sweeps[cell_idx_path[i+1]])){
                    next_point = findNextGoal(point, cells_sweeps[cell_idx_path[i+1]].back(), next_candidates);
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]], point, next_point);
                    way_points.insert(way_points.end(), std::next(shortest_path.begin()), std::prev(shortest_path.end()));
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i+1]], next_point, cells_sweeps[cell_idx_path[i+1]].back());
                }else{
                    next_point = findNextGoal(point, cells_sweeps[cell_idx_path[i+1]].front(), next_candidates);
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]], point, next_point);
                    way_points.insert(way_points.end(), std::next(shortest_path.begin()), std::prev(shortest_path.end()));
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i+1]], next_point, cells_sweeps[cell_idx_path[i+1]].front());
                }
                way_points.insert(way_points.end(), shortest_path.begin(), std::prev(shortest_path.end()));
                point = way_points.back();
            }
        }else{
            shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]],
                                            cells_sweeps[cell_idx_path[i]].front(),
                                            cells_sweeps[cell_idx_path[i]].back());
            if(doReverseNextSweep(point, cells_sweeps[cell_idx_path[i]])){
                way_points.insert(way_points.end(), shortest_path.rbegin(), shortest_path.rend());
            }else{
                way_points.insert(way_points.end(), shortest_path.begin(), shortest_path.end());
            }
            point = way_points.back();

            if((i+1)<cell_idx_path.size()){
                next_candidates = cell_intersections[cell_idx_path[i]][cell_idx_path[i+1]];
                if(doReverseNextSweep(point, cells_sweeps[cell_idx_path[i+1]])){
                    next_point = findNextGoal(point, cells_sweeps[cell_idx_path[i+1]].back(), next_candidates);
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]], point, next_point);
                    way_points.insert(way_points.end(), std::next(shortest_path.begin()), std::prev(shortest_path.end()));
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i+1]], next_point, cells_sweeps[cell_idx_path[i+1]].back());
                }else{
                    next_point = findNextGoal(point, cells_sweeps[cell_idx_path[i+1]].front(), next_candidates);
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]], point, next_point);
                    way_points.insert(way_points.end(), std::next(shortest_path.begin()), std::prev(shortest_path.end()));
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i+1]], next_point, cells_sweeps[cell_idx_path[i+1]].front());
                }
                way_points.insert(way_points.end(), shortest_path.begin(), std::prev(shortest_path.end()));
                point = way_points.back();
            }
        }
    }

    cv::Mat img(600, 600, CV_8UC3, cv::Scalar(255, 255, 255));
    // draw poly contours
    cv::Point p1, p2;
    cv::namedWindow("cover",cv::WINDOW_NORMAL);
    cv::imshow("cover", img);

    for (int i = 1; i < ext_poly.size(); i++)
    {
        p1 = cv::Point(CGAL::to_double(ext_poly[i-1].x()),CGAL::to_double(ext_poly[i-1].y()));
        p2 = cv::Point(CGAL::to_double(ext_poly[i].x()),CGAL::to_double(ext_poly[i].y()));
        cv::line(img, p1, p2, cv::Scalar(0, 0, 0));
        cv::namedWindow("cover",cv::WINDOW_NORMAL);
        cv::imshow("cover", img);
        cv::line(img, p1, p2, cv::Scalar(0, 0, 0));
    }
    for (int i = 0; i < inner_polys.size(); i++)
    {
        inner_polys[i].emplace_back(inner_polys[i].front());
        for (int j = 1; j < inner_polys[i].size(); j++)
        {
            p1 = cv::Point(CGAL::to_double(inner_polys[i][j-1].x()),CGAL::to_double(inner_polys[i][j-1].y()));
            p2 = cv::Point(CGAL::to_double(inner_polys[i][j].x()),CGAL::to_double(inner_polys[i][j].y()));
            cv::line(img, p1, p2, cv::Scalar(255, 0, 0));
            cv::namedWindow("cover",cv::WINDOW_NORMAL);
            cv::imshow("cover", img);
            cv::line(img, p1, p2, cv::Scalar(255, 0, 0));
        }
    }
    // draw waypoints
    for(size_t i = 1; i < way_points.size(); ++i){
        p1 = cv::Point(CGAL::to_double(way_points[i-1].x()),CGAL::to_double(way_points[i-1].y()));
        p2 = cv::Point(CGAL::to_double(way_points[i].x()),CGAL::to_double(way_points[i].y()));
        cv::line(img, p1, p2, cv::Scalar(0, 64, 255));
        cv::namedWindow("cover",cv::WINDOW_NORMAL);
        cv::imshow("cover", img);
        cv::waitKey(50);
        cv::line(img, p1, p2, cv::Scalar(200, 200, 200));
    }

    std::cout << "Output waypoint size: " << way_points.size() << std::endl;
    return way_points;
}

void contourCallback(const nav_msgs::Path::ConstPtr& path_msg)
{
    if (path_msg->poses.size() != 0)
    {
        std::vector<Point_2> poly_contour;
        for (int i = 0; i < path_msg->poses.size(); i++)
        {
            //std::cout << path_msg->poses[i].pose.position.x << ", " << path_msg->poses[i].pose.position.y << std::endl;
            poly_contour.emplace_back(Point_2(path_msg->poses[i].pose.position.x, path_msg->poses[i].pose.position.y));
        }
        polys.emplace_back(poly_contour);
        ROS_INFO("Received contour %u with %zu poses.", contour_rcv_count, path_msg->poses.size());
        contour_rcv_count++;
    }
    else
    {
        ROS_INFO("Start Planning global path");
        globalPathPlan(polys, start);
    }
}

int main(int argc, char **argv)
{
    // 初始化ROS节点
    ros::init(argc, argv, "global_planner");
    ros::NodeHandle nh;

    ros::Subscriber contour_sub = nh.subscribe("contour", 10, contourCallback);

    ros::spin();

    return 0;
}
