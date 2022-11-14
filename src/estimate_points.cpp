#include <Eigen/Dense>

#include <iostream>
#include <fstream>

#include <string>

#include <vector>
#include <map>

#include <algorithm>

#include <math.h>


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
    A << slope_i,1,
         slope_j,1;

    Eigen::Matrix<double,2,1> b;
    b << intercept_i,
         intercept_j;

    Eigen::Matrix<double,2,1> sol;

    sol = A.inverse() * b;

    x = sol(0,0);
    y = sol(1,0);
}


// find intersection between all the points and take 
void estimate_points_greedy(std::vector<std::vector<double>> lines, int n_points)
{
    // coordinates of intersection and number of intersections there
    std::map<std::pair<double,double>,int> points;


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

                std::pair point = std::make_pair(x_int,y_int);

                if(points.find(point) == points.end())
                {
                    points[point] = 1;
                }
                else
                {
                    ++points[point];

                    std::cout << "point.first: " << point.first << " point.second: " << point.second << "| n: " << points[point] << std::endl;
                }

            }

        }

    }

    auto pr = std::max_element( points.begin(), points.end(), [](const auto &x, const auto &y)
                            {
                                return x.second < y.second;
                            });

    std::cout << "pr->first: " << pr->first.first << ", " << pr->first.second << std::endl;
    std::cout << "pr->second: " << pr->second << std::endl;

    std::cout << "n_intersections: " << lines.size() * lines.size() << std::endl;
    std::cout << "points.size(): " << points.size() << std::endl;


}

int main()
{

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
    inFile.open("../Task-Points/1-lines.txt");
    
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

    estimate_points_greedy(lines,n_points);

    std::cout << "n_lines: " << n_lines << "| n_points: " << n_points << std::endl;

    std::cout << "lines.size(): " << lines.size() << std::endl;




    return 0;
}