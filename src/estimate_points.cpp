#include <Eigen/Dense>

#include <iostream>
#include <fstream>

#include <string>

#include <vector>
#include <map>

#include <algorithm>

#include <math.h>

#include <chrono>

bool cmp(std::pair<std::pair<double,double>,int> &a,
         std::pair<std::pair<double,double>,int> &b)
{
    return a.second < b.second;
}

void sort(std::map<std::pair<double,double>,int> &my_map, std::vector<std::pair<std::pair<double,double>,int>> &sorted_map)
{
    // std::vector<std::pair<std::pair<double,double>,int>> sorted_map;

    for(auto &it : my_map)
    {
        sorted_map.push_back(it);
    }

    std::sort(sorted_map.begin(), sorted_map.end(), cmp);

    // for(auto &it: sorted_map)
    // {
    //     std::cout << it.first.first << " " << it.first.second << " | " << it.second << std::endl;
    // }
}


double euclidean_distance(const std::pair<double, double> &pt1, const std::pair<double,double> &pt2)
{
    return sqrt((pt1.first-pt2.first)*(pt1.first-pt2.first) + (pt1.second-pt2.second)*(pt1.second-pt2.second));
}

// in task there is max allowed displacement from cluster 
void find_closest_cluster(std::vector<std::pair<double,double>> clusters, std::pair<double,double> point, int &cluster_num, double max_dist = 5)
{
    for(size_t i=0; i<clusters.size(); i++)
    {
        double dist = euclidean_distance(point, clusters[i]);

        if(dist < max_dist)
        {
            cluster_num = i;
            return;
        }

    }

    cluster_num = -1;
}


// clustering intersection points
// with known number of points we only have to find N points with highest intersection number
// input is already sorted
void cluster_points(std::vector<std::pair<std::pair<double,double>,int>> sorted_map, std::vector<std::pair<double, double>> &clusters, std::vector<int> &cnts, int n_clusters)
{
    int n_clusters_found = 0;


    // std::vector<std::pair<double,double>> clusters;
    // std::vector<int> cnts;

    // clusters.resize(n_clusters);
    // cnts.resize(n_clusters);

    for(size_t i = sorted_map.size()-1; i>=0; --i)
    {
        if(n_clusters==n_clusters_found)
        {
            break;
        }

        // initialize first cluster
        if(n_clusters_found==0)
        {
            std::pair<double, double> pt(sorted_map[i].first);


            clusters.push_back(pt);
            cnts.push_back(sorted_map[i].second);
            // clusters[n_clusters_found] = pt;
            // cnts[n_clusters_found] = sorted_map[i][1];
            
            n_clusters_found++;
        }

        else
        {
            // get coordinates of the point
            std::pair<double, double> pt(sorted_map[i].first);
            // find if it is already close enough to some cluster:
            int cluster_num = -1;

            find_closest_cluster(clusters, pt, cluster_num);

            // if could not find point that is closer than threshold -> create new cluster
            if(cluster_num==-1)
            {
                clusters.push_back(pt);
                cnts.push_back(sorted_map[i].second);
                ++n_clusters_found;
            }
            else
            {
                cnts[cluster_num] += sorted_map[i].second;
            }
            // if close point was found -> point belongs to that cluster
        }



    }


}

std::vector<std::string> split_by_char(std::string input_line, char splitter)
{
    std::vector<std::string> splitted_input;

    std::string input_cur;
    input_cur.clear();
    for(auto x : input_line)
    {
        // std::cout << "x: " << x << "input_cur: " << input_cur << std::endl;
        if(x == splitter)
        {
            splitted_input.push_back(input_cur);
            input_cur.clear();
        }
        else
        {
            input_cur += x;
        }
    }

    splitted_input.push_back(input_cur);

    return splitted_input;
}

// estimate intercept and slope of a line y = mx+b,
// b - intercept, m - slope
void get_line_coefs(std::vector<double> line, double &intercept, double &slope)
{
    double x_1,y_1,x_2,y_2;

    x_1 = line[0];
    y_1 = line[1];
    x_2 = line[2];
    y_2 = line[3];

    slope = (y_2 - y_1)/(x_2 - x_1);
    intercept = y_1 - slope*x_1;

}


// estimate point of intersection between 2 lines
void get_lines_intersection(double intercept_i, double slope_i, double intercept_j, double slope_j, double &x, double &y)
{
    // solve via matrix inverse
    Eigen::Matrix<double,2,2> A;
    A << slope_i,-1,
         slope_j,-1;

    Eigen::Matrix<double,2,1> b;
    b << -intercept_i,
         -intercept_j;

    Eigen::Matrix<double,2,1> sol;

    sol = A.inverse() * b;

    x = sol(0,0);
    y = sol(1,0);
}

void get_lines_intersection_simple(double intercept_i, double slope_i, double intercept_j, double slope_j, double &x, double &y)
{
    x = (intercept_i - intercept_j)/(slope_j-slope_i);
    y = slope_j*x + intercept_j;
}


// find intersection between all the points and take 
void estimate_points_greedy(std::vector<std::vector<double>> lines, int n_points)
{
    // coordinates of intersection and number of intersections there
    std::map<std::pair<double,double>,int> points;

    line_coefs

    auto intersection_start = std::chrono::high_resolution_clock::now();

    for(size_t i = 0; i < lines.size(); i++)
    {
        double intercept_i, slope_i;
        get_line_coefs(lines[i], intercept_i, slope_i);
        for(size_t j = 0; j < lines.size(); j++)
        {
            if(i != j)
            {
                double intercept_j, slope_j;
                get_line_coefs(lines[j], intercept_j, slope_j);

                double x,y;
                get_lines_intersection(intercept_i,slope_i, intercept_j, slope_j, x,y);

                // std::cout << "x: " << x << "y: " << y << std::endl;


                int x_int = round(x);
                int y_int = round(y);


                if(x>=-1000 && x<= 1000 && y>=-1000 && y<=1000)
                {
                                    // intersection_duration: 866437 milliseconds (original)

                    std::pair point = std::make_pair(x_int,y_int);

                    if(points.find(point) == points.end())
                    {
                        points[point] = 1;
                    }
                    else
                    {
                        ++points[point];
                        // std::cout << "point.first: " << point.first << " point.second: " << point.second << "| n: " << points[point] << std::endl;
                    }
                }



            }

        }

    }

    auto intersection_end = std::chrono::high_resolution_clock::now();
    auto intersection_duration = std::chrono::duration_cast<std::chrono::milliseconds>(intersection_end - intersection_start);
    std::cout << "intersection_duration: " << intersection_duration.count() << " milliseconds" << std::endl;


    std::vector<std::pair<std::pair<double,double>,int>> sorted_map;


    auto sort_start = std::chrono::high_resolution_clock::now();
    sort(points, sorted_map);

    auto sort_end = std::chrono::high_resolution_clock::now();

    auto sort_duration = std::chrono::duration_cast<std::chrono::milliseconds>(sort_end - sort_start);
    std::cout << "sort_duration: " << sort_duration.count() << " milliseconds" << std::endl;



    auto cluster_start = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<double,double>> clusters;
    std::vector<int> cnts;
    cluster_points(sorted_map, clusters, cnts, n_points);

    auto cluster_end = std::chrono::high_resolution_clock::now();
    auto cluster_duration = std::chrono::duration_cast<std::chrono::milliseconds>(cluster_end - cluster_start);
    std::cout << "cluster_duration: " << cluster_duration.count() << " milliseconds" << std::endl;

    for(size_t i=0;i<clusters.size();i++)
    {
        std::cout << "pt.x: " << clusters[i].first << ", pt.y: " << clusters[i].second << " count: " << cnts[i] << std::endl;
    }


    // auto pr = std::max_element( points.begin(), points.end(), [](const auto &x, const auto &y)
    //                         {
    //                             return x.second < y.second;
    //                         });

    // std::cout << "pr->first: " << pr->first.first << ", " << pr->first.second << std::endl;
    // std::cout << "pr->second: " << pr->second << std::endl;

    // std::cout << "n_intersections: " << lines.size() * lines.size() << std::endl;
    // std::cout << "points.size(): " << points.size() << std::endl;


}

int main()
{
    double a1 = 2;
    double b1 = 4;
    double a2 = -2;
    double b2 = 2;

    double x_point,y_point;
    
    get_lines_intersection_simple(b1,a1,b2,a2,x_point,y_point);

    std::cout<< "INTERSECTION. x: " << x_point << " y: " << y_point << std::endl;

    Eigen::Matrix<double,2,2> A;
    A(0,0) = 1.0;
    A(0,1) = 2.0;
    A(1,0) = 2.0;
    A(1,1) = 3.0;

    Eigen::Matrix<double,2,1> b;
    b(0,0) = -1.0;
    b(1,0) = -5.0;


    Eigen::Matrix<double,2,1> x;

    x = A.inverse() * b;

    std::cout << x << std::endl;


    std::ifstream inFile;
    std::string inLine;
    inFile.open("../Task-Points/3-lines.txt");
    
    bool first_line = true;
    int n_points = 0;
    int n_lines = 0;
    
    std::vector<std::string> input_splitted;
    // lines as 2 points (x_1,y_1), (x_2,y_2)
    std::vector<std::vector<double>> lines;
    while(std::getline(inFile, inLine))
    {
        // split input by space
        input_splitted = split_by_char(inLine,' ');
        if(first_line)
        {
            first_line = false;
            n_lines = std::stod(input_splitted[0]);
            n_points = std::stod(input_splitted[1]);

        }

        

        else
        {
            double x_1 = std::stod(input_splitted[0]);
            double y_1 = std::stod(input_splitted[1]);
            double x_2 = std::stod(input_splitted[2]);
            double y_2 = std::stod(input_splitted[3]);

            std::vector<double> line = {x_1,y_1,x_2,y_2};

            lines.push_back(line);

        }

    }

    auto estimation_start = std::chrono::high_resolution_clock::now();
    estimate_points_greedy(lines,n_points);
    auto estimation_end = std::chrono::high_resolution_clock::now();

    auto estimation_duration = std::chrono::duration_cast<std::chrono::milliseconds>(estimation_end - estimation_start);

    std::cout << "estimation_duration: " << estimation_duration.count() << " milliseconds" << std::endl;

    std::cout << "n_lines: " << n_lines << "| n_points: " << n_points << std::endl;

    std::cout << "lines.size(): " << lines.size() << std::endl;




    return 0;
}