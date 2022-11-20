
#include <algorithm>
#include <math.h>


#include <omp.h>
#define THREAD_NUM 4

#include <iomanip>
#include <iostream>

#include <vector>
#include <map>

#include <chrono>
#include <fstream>
// #include <Eigen/Dense>


// simple euclidean norm between 2 points
double euclidean_distance(const std::pair<double, double> &pt1, const std::pair<double,double> &pt2)
{
    return sqrt((pt1.first-pt2.first)*(pt1.first-pt2.first) + (pt1.second-pt2.second)*(pt1.second-pt2.second));
}


// compare function to sort map by values
bool cmp(std::pair<std::pair<double,double>,int> &a,
         std::pair<std::pair<double,double>,int> &b)
{
    return a.second < b.second;
}

/*
    sorts map by values and writes them into vector
*/
void sort(std::map<std::pair<double,double>,int> &my_map, std::vector<std::pair<std::pair<double,double>,int>> &sorted_map)
{

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



/*
    checking if point belongs to existing cluster
    if does not -> new cluster will be created
*/ 
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


/*
    clustering intersection points
    with known number of points we only have to find N points with highest intersection number
    input is already sorted
*/ 
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

/*
    estimate intercept and slope of a line y = mx+b,
    b - intercept, m - slope
*/ 
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

// estimate point of intersection between 2 lines with linear algebra
// void get_lines_intersection(double intercept_i, double slope_i, double intercept_j, double slope_j, double &x, double &y)
// {
//     // solve via matrix inverse
//     Eigen::Matrix<double,2,2> A;
//     A << slope_i,-1,
//          slope_j,-1;

//     Eigen::Matrix<double,2,1> b;
//     b << -intercept_i,
//          -intercept_j;

//     Eigen::Matrix<double,2,1> sol;

//     sol = A.inverse() * b;

//     x = sol(0,0);
//     y = sol(1,0);
// }

// estimate point of intersection by simple substitution(works faster)
void get_lines_intersection_simple(double intercept_i, double slope_i, double intercept_j, double slope_j, double &x, double &y)
{
    x = (intercept_i - intercept_j)/(slope_j-slope_i);
    y = slope_j*x + intercept_j;
}


/*
    finds all points from line intersections
    1. calculate all line coeffs(slope, intercept)
    2. calculate intersection points(rounding to closest integer) and number of intersections
    3. sort intersection points by number of intersections
    4. cluster intersection points(otherwise there are several points very close to real point)

*/ 
void estimate_points_greedy(std::vector<std::vector<double>> lines, int n_points, bool print_clusters=false)
{
    // coordinates of intersection and number of intersections there
    std::map<std::pair<double,double>,int> points;

    
    auto intersection_start = std::chrono::high_resolution_clock::now();


    // precalculate all lines coeffs(intercept, slope)
    std::vector<std::pair<double,double>> line_coefs;
    line_coefs.resize(lines.size());
    double intercept, slope;
    for(size_t i=0; i<lines.size(); i++)
    {
        get_line_coefs(lines[i], intercept, slope);
        line_coefs[i] = std::pair(intercept,slope);
    }

    omp_set_num_threads(THREAD_NUM);
    #pragma omp parallel for
    for(size_t i = 0; i < lines.size(); i++)
    {
        double intercept_i, slope_i;
        // get_line_coefs(lines[i], intercept_i, slope_i);
        intercept_i = line_coefs[i].first;
        slope_i = line_coefs[i].second;
        for(size_t j = 0; j < lines.size(); j++)
        {
            if(i != j)
            {
                double intercept_j, slope_j;
                intercept_j = line_coefs[j].first;
                slope_j = line_coefs[j].second;
                // get_line_coefs(lines[j], intercept_j, slope_j);

                double x,y;
                get_lines_intersection_simple(intercept_i, slope_i, intercept_j, slope_j, x,y);

                // std::cout << "x: " << x << "y: " << y << std::endl;


                int x_int = round(x);
                int y_int = round(y);


                if(x>=-1000 && x<= 1000 && y>=-1000 && y<=1000)
                {
                    // intersection_duration: 866437 milliseconds (original)
                    // estimation_duration: 517210 milliseconds(after reducing range to [-1000,1000])
                    // duration after using simple intersection calculation instead of linalg: 354135 milliseconds
                    // duration after precomputing line coeffs: 268880 milliseconds
                    // after adding multithreading into 4 threads: 102222

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

    if(print_clusters)
    {
        for(size_t i=0;i<clusters.size();i++)
        {
            std::cout << "pt.x: " << clusters[i].first << ", pt.y: " << clusters[i].second << " count: " << cnts[i] << std::endl;
        }
    }



    // writing clusters to file:
    std::ofstream outFile;
    outFile.open("out.txt");
    
    outFile << std::setprecision(2) << std::fixed;
    for(size_t i=0; i<clusters.size(); i++)
    {
        
        outFile << static_cast<double>(clusters[i].first) << " " << clusters[i].second << "\n"; 
    }
    outFile.close();



}